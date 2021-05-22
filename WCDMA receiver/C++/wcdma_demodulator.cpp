#include "lime/LimeSuite.h"
#include <iostream>
#include<vector>
#include <memory>
#include<complex>
#include <fstream>
#include<conio.h>
#include "ipps.h"

using namespace std;

//Device structure, should be initialize to NULL
lms_device_t* dev1 = NULL;

int error()
{
	if (dev1 != NULL)
		LMS_Close(dev1);
	exit(-1);
}

class wcdma_receiver
{
	typedef std::complex<float> IQf;
public:
	wcdma_receiver(int sps_log2, float f_d, int decimate_rate, int poly_bank_size,
		float dll_bw, float damp_ratio,ofstream &out) :f_d{ f_d }, decimate_rate{ decimate_rate }, 
		sps_log2{ sps_log2 }, poly_bank_size{ poly_bank_size }, out{out}
	{
		sps = 1 << sps_log2;
		FFT_SIZE = 1 << (sps_log2 + 16);
		samples_in_frame = sps * NCHIPS_IN_FRAME;
		int SpecSize, SpecBufferSize, BufferSize;
		IppStatus status = ippsFFTGetSize_C_32fc(sps_log2 + 16,
			IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, &SpecSize, &SpecBufferSize, &BufferSize);

		SpecBuffer.resize(SpecSize), init_space.resize(SpecBufferSize), wspace.resize(BufferSize);

		status = ippsFFTInit_C_32fc(&pFFTSpec, sps_log2 + 16,
			IPP_FFT_DIV_INV_BY_N, ippAlgHintNone, SpecBuffer.data(), init_space.data());

		pattern.resize(FFT_SIZE), sig_fc.resize(FFT_SIZE), mags.resize(FFT_SIZE),
			frame_fc.resize(FFT_SIZE), frame_chips.resize(FFT_SIZE);
		circle_buf.resize(BUF_SZ), WHT.resize(10);

		for (int i = 0; i < 10; i++)
			WHT[i].resize(NCHIPS_IN_FRAME);

		ir_sz = RRC_HALF_IR_LEN_IN_CHIPS * 2 * sps + 1;
		float omega_n = 8 * damp_ratio * dll_bw / (4 * damp_ratio * damp_ratio + 1);
		const float Td = CPICH_SF * sps / f_d;
		float A = omega_n * Td;
		k0 = (8 * damp_ratio * A) / (4 + 4 * damp_ratio * A + A * A);
		k1 = 4 * A * A / (4 + 4 * damp_ratio * A + A * A);

		bank.resize(poly_bank_size);
		for (int i = 0; i < poly_bank_size; i++)
			bank[i].resize(ir_sz);

		for (int i = 0; i < poly_bank_size; i++)
			for (int j = 0; j < ir_sz; j++)
				bank[i][ir_sz - 1 - j] = { rrc(alpha, float(i + j * poly_bank_size) / decimate_rate -
					RRC_HALF_IR_LEN_IN_CHIPS),0 };

		state = ESTIMATIONS;
	}

	int general_work(std::vector<Ipp16sc> &sig,bool &slow)
	{
		int psc_code = 0, start_pos;
		std::vector<IQf> cpich_;
		float m;
		switch (state)
		{
		case ESTIMATIONS:
			freq_correction_val = measure_initial_df(sig);
			psc_code = estimate_primary_code(sig, freq_correction_val);

			std::cout << "psc code: " << psc_code << std::endl;

			cpich = gen_cpich(16 * psc_code);

			diff_cpich = {};
			for (int i = 0; i < cpich.size(); i++)
			{
				IQf elem = cpich[(i + 1) % NCHIPS_IN_FRAME] -
					cpich[(i - 1 + NCHIPS_IN_FRAME) % NCHIPS_IN_FRAME];
				diff_cpich.push_back(std::conj(elem));
			}

			cpich_ = upsample(cpich, sps);
			form_pattern_in_freq_domain(cpich_);

			ippsConj_32fc_I((Ipp32fc*)cpich.data(), NCHIPS_IN_FRAME);

			state = CATCH;
			return FFT_SIZE;

		case CATCH:
			ippsConvert_16s32f((Ipp16s*)sig.data(), (Ipp32f*)sig_fc.data(), FFT_SIZE * 2);
			establish_initial_state(freq_correction_val);
			res(sig_fc, false);
			ippsCopy_32fc((Ipp32fc*)frame_chips.data(), sig_fc.data(), FFT_SIZE);
			state = PASS_SAMPLES;
			return samples_in_frame - ((FFT_SIZE - fast_conv(sig_fc, m)) % samples_in_frame);

		case PASS_SAMPLES:
			state = WORK;
			establish_initial_state(freq_correction_val);
			return FFT_SIZE;

		case WORK:
			ippsConvert_16s32f((Ipp16s*)sig.data(), (Ipp32f*)sig_fc.data(), FFT_SIZE * 2);
			res(sig_fc,true);
			if (state == ESTIMATIONS)
				slow = true;
			return FFT_SIZE;
		}
	}

private:
	std::vector<Ipp32fc>pattern, sig_fc, frame_fc;
	std::vector<Ipp32f>mags;
	std::vector<Ipp8u> SpecBuffer, init_space, wspace, firbuf;
	IppsFFTSpec_C_32fc *pFFTSpec;
	int sps,sps_log2, FFT_SIZE, samples_in_frame,nlosts_counter;
	float f_d, freq_correction_val, k0, k1, dll_loop, remain_frac_time_error;
	enum { ESTIMATIONS, CATCH, PASS_SAMPLES, WORK } state;
	const int CYCLE_LEN = 262143;
	const int NCHIPS_IN_FRAME = 38400;
	const int RRC_HALF_IR_LEN_IN_CHIPS = 2;
	const int BUF_SZ = 1000000;
	const int CPICH_SF_LOG2 = 8;
	const int CPICH_SF = 1 << CPICH_SF_LOG2;
	const int cpich_errors_threshold = 100;
	const float PI = 3.1415926;
	const float Tframe = 0.01f;
	const float Kd = 1.9;
	const float df_max = 2500;
	const float df_step = 20;
	const float alpha = 0.22f;
	std::vector<IQf> cpich, frame_chips, diff_cpich, circle_buf;
	uint64_t free, last_sample_used_in_res,packet_num;
	IQf phase, theta;
	std::vector<std::vector<IQf>> WHT;
	std::vector<std::vector<Ipp32fc>> bank;
	int poly_bank_size, decimate_rate, ir_sz, bank_id, id_of_chip_to_get,
		how_many_samples_need_to_read;
	Ipp32fc previous_weight_estimation;
	ofstream &out;

	int fast_conv(std::vector<Ipp32fc> &sig_fc, float &m)
	{
		int start_pos;
		ippsFFTFwd_CToC_32fc_I(sig_fc.data(), pFFTSpec, wspace.data());
		ippsMul_32fc_I(pattern.data(), sig_fc.data(), FFT_SIZE);
		ippsFFTInv_CToC_32fc_I(sig_fc.data(), pFFTSpec, wspace.data());
		ippsMagnitude_32fc(sig_fc.data(), mags.data(), FFT_SIZE);
		ippsMaxIndx_32f(mags.data(), FFT_SIZE, &m, &start_pos);
		return start_pos;
	}

	void form_pattern_in_freq_domain(std::vector<IQf> &time_domain)
	{
		ippsZero_32fc(pattern.data(), FFT_SIZE);
		ippsCopy_32fc((Ipp32fc*)time_domain.data(), pattern.data(), time_domain.size());
		ippsFFTFwd_CToC_32fc_I(pattern.data(), pFFTSpec, wspace.data());
		ippsConj_32fc_I(pattern.data(), FFT_SIZE);
	}

	int estimate_primary_code(std::vector<Ipp16sc> &sig, float df)
	{
		float magnitude = -1, m;
		int psc_code, start_pos;

		ippsConvert_16s32f((Ipp16s*)sig.data(), (Ipp32f*)sig_fc.data(), FFT_SIZE * 2);
		establish_initial_state(df);
		res(sig_fc, false);
		ippsCopy_32fc((Ipp32fc*)frame_chips.data(), frame_fc.data(), FFT_SIZE);

		for (int i = 0; i < 512; i++)
		{
			ippsCopy_32fc(frame_fc.data(), sig_fc.data(), FFT_SIZE);

			std::vector<IQf> cpich = gen_cpich(16 * i);
			cpich = upsample(cpich, sps);
			form_pattern_in_freq_domain(cpich);

			fast_conv(sig_fc, m);

			if (m > magnitude)
			{
				magnitude = m;
				psc_code = i;
			}
		}

		return psc_code;
	}

	float measure_initial_df(std::vector<Ipp16sc> &sig)
	{
		float sig_df = 0, magnitude = -1, m;

		std::vector<IQf> psc = upsample(generate_psc_pattern(), sps);
		form_pattern_in_freq_domain(psc);

		for (float df = -df_max; df < df_max; df += df_step)
		{
			ippsConvert_16s32f((Ipp16s*)sig.data(), (Ipp32f*)sig_fc.data(), FFT_SIZE * 2);			
			establish_initial_state(df);
			res(sig_fc,false);
			ippsCopy_32fc((Ipp32fc*) frame_chips.data(), sig_fc.data(),FFT_SIZE);
			
			fast_conv(sig_fc, m);

			if (m > magnitude)
			{
				magnitude = m;
				sig_df = df;
			}
		}

		return sig_df;
	}

	float rrc(float a, float t)
	{
		t += 0.000001f;
		return (sin(PI * t * (1 - a)) + 4 * a * t * cos(PI * t * (1 + a)))
			/ (PI * t * (1 - 16 * a * a * t * t));
	}

	std::vector<IQf> gen_cpich(int p)
	{
		std::vector<float> L1 = { -1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
		std::vector<float> L2 = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };
		std::vector<float> M1(CYCLE_LEN), M2(CYCLE_LEN), GOLD(CYCLE_LEN);
		std::vector<IQf> CPICH(NCHIPS_IN_FRAME);

		for (int i = 0; i < CYCLE_LEN; i++)
		{
			M1[i] = L1[0], M2[i] = L2[0];
			float m1 = L1[0] * L1[7];
			float m2 = L2[0] * L2[5] * L2[7] * L2[10];

			for (int j = 0; j < 17; j++)
				L1[j] = L1[j + 1], L2[j] = L2[j + 1];
			L1[17] = m1, L2[17] = m2;
		}

		for (int i = 0; i < CYCLE_LEN; i++)
			GOLD[i] = M1[(i + p) % CYCLE_LEN] * M2[i];

		for (int i = 0; i < NCHIPS_IN_FRAME; i++)
			CPICH[i] = { GOLD[i],GOLD[(i + CYCLE_LEN / 2 + 1) % CYCLE_LEN] };

		return CPICH;
	}

	std::vector<IQf> upsample(std::vector<IQf> in, int sps)
	{
		std::vector<IQf> out;
		for (auto e : in)
		{
			out.push_back(e);
			for (int i = 1; i < sps; i++)
				out.push_back(0);
		}

		return out;
	}

	std::vector<IQf> generate_psc_pattern()
	{
		std::vector<float> a = { 1, 1, 1, 1, 1, 1, -1, -1, 1, -1, 1,-1,1,-1,-1, 1 };
		std::vector<float> cp = { 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1 };

		std::vector<IQf> psc;
		for (int slot = 0; slot < 15; slot++)
		{
			for (int m = 0; m < 16; m++)
				for (int n = 0; n < 16; n++)
					psc.push_back({ cp[m] * a[n], cp[m] * a[n] });
			for (int chip = 0; chip < 2304; chip++)
				psc.push_back({ 0,0 });
		}

		return psc;
	}

	void establish_initial_state(float freq_correction_val)
	{
		packet_num = 0;
		dll_loop = remain_frac_time_error = 0;
		phase = 1;
		theta = std::polar(1.f, freq_correction_val * 2 * PI / f_d);
		free = last_sample_used_in_res = 0;
		id_of_chip_to_get = how_many_samples_need_to_read = 0;
		bank_id = 0;
		previous_weight_estimation = { 0,0 };

		for (int i = 0; i < BUF_SZ; i++)
			circle_buf[i] = 0;
	}

	void res(std::vector<Ipp32fc> &sig, bool mode)
	{
		for (int i = 0; i < FFT_SIZE; i++)
		{
			circle_buf[free % BUF_SZ] = { sig[i].re, sig[i].im };
			free++;
		}

		do
		{
			for (int i = 0; i < how_many_samples_need_to_read; i++)
			{
				last_sample_used_in_res++;
				circle_buf[last_sample_used_in_res % BUF_SZ] *= phase;
				phase *= theta;
			}

			frame_chips[id_of_chip_to_get] = 0;
			int pos = (last_sample_used_in_res + BUF_SZ - ir_sz + 1) % BUF_SZ;
			if (pos + ir_sz - 1< BUF_SZ)
				ippsDotProd_32fc(bank[bank_id].data(), (Ipp32fc*)circle_buf.data() + pos,ir_sz, 
				(Ipp32fc*)frame_chips.data() + id_of_chip_to_get);

			if (id_of_chip_to_get % CPICH_SF == CPICH_SF - 1 && mode)
			{
				int shift = id_of_chip_to_get - CPICH_SF + 1;
				Ipp32fc time_err, ampl;
				ippsDotProd_32fc((Ipp32fc*)frame_chips.data() + shift, (Ipp32fc*)diff_cpich.data() + shift, CPICH_SF, &time_err);
				ippsDotProd_32fc((Ipp32fc*)frame_chips.data() + shift, (Ipp32fc*)cpich.data() + shift, CPICH_SF, &ampl);

				if (ampl.re != 0 || ampl.im != 0)
				{
					float dt = -std::real(IQf{ time_err.re,time_err.im } / IQf{ ampl.re,ampl.im }) / Kd;
					dll_loop += dt * k1;
					float control_signal = k0 * dt + dll_loop;
					remain_frac_time_error += control_signal * decimate_rate;
					bank_id += int(remain_frac_time_error);
					remain_frac_time_error = remain_frac_time_error - int(remain_frac_time_error);
				}
			}
			id_of_chip_to_get++;

			if (mode)
				bank_id += decimate_rate;
			else
				bank_id += decimate_rate / sps;
			how_many_samples_need_to_read = bank_id / poly_bank_size;
			bank_id %= poly_bank_size;

			if (id_of_chip_to_get == NCHIPS_IN_FRAME && mode) //start to process a frame
			{

				id_of_chip_to_get = 0;

				Ipp32fc weight;
				ippsMul_32fc_I((Ipp32fc*)cpich.data(), (Ipp32fc*)frame_chips.data(), NCHIPS_IN_FRAME);
				ippsSum_32fc((Ipp32fc*)frame_chips.data(), NCHIPS_IN_FRAME, &weight, ippAlgHintFast);

				// frequency error estimation:
				float re = previous_weight_estimation.re;
				float im = previous_weight_estimation.im;
				float error_angle = 0;
				if ((re != 0 || im != 0) && (weight.re != 0 && weight.im != 0))
					error_angle = arg(IQf{ weight.re,weight.im } *IQf{ re,-im });
				freq_correction_val -= error_angle * 0.01f / (Tframe * PI * 2);
				theta = std::polar(1.f, freq_correction_val * 2 * PI / f_d);
				phase /= abs(phase);
				previous_weight_estimation = weight;

				ippsMulC_32fc_I({ weight.re,-weight.im }, (Ipp32fc*)frame_chips.data(), NCHIPS_IN_FRAME);
				ippsMulC_32fc_I({ 1,1 }, (Ipp32fc*)frame_chips.data(), NCHIPS_IN_FRAME);

				// walsh-hadamard
				ippsCopy_32fc((Ipp32fc*)frame_chips.data(), (Ipp32fc*)WHT[0].data(), NCHIPS_IN_FRAME);
				for (int i = 1; i < 10; i++)
					for (int j = 0; j < NCHIPS_IN_FRAME/2; j++)
					{
						WHT[i][j] = WHT[i - 1][j * 2] + WHT[i - 1][j * 2 + 1];
						WHT[i][j + NCHIPS_IN_FRAME/2] = WHT[i - 1][j * 2] - WHT[i - 1][j * 2 + 1];
					}
				
				// control synch mode
				int errors_in_cpich = 0;
				for (int i = 0; i < NCHIPS_IN_FRAME / CPICH_SF; i++)
				{
					if (WHT[8][i].real() < 0) errors_in_cpich++;
					if (WHT[8][i].imag() < 0) errors_in_cpich++;
				}

				if (errors_in_cpich > cpich_errors_threshold)
				{
					std::cout << "SIGNAL WAS LOST!!!" << std::endl;
					
					state = ESTIMATIONS;
				}

				if (packet_num++ % 10 == 0)
				{
					std::cout << "errors in cpich: " << errors_in_cpich << std::endl;
					//CPICH SNR estimation:
					Ipp32f A, ro;
					ippsMeanStdDev_32f((Ipp32f*)WHT[8].data(), NCHIPS_IN_FRAME / CPICH_SF * 2, &A, &ro, ippAlgHintFast);
					if (ro != 0)
						std::cout << "input CPICH SNR: " << 20 * log10(A / ro) - 3.01 * CPICH_SF_LOG2
						- 3.01 * sps_log2 - 3.01 << std::endl;

					std::cout << "freq error: " << -freq_correction_val << std::endl;
					//out << -freq_correction_val << std::endl;
				}
			}
		} while (last_sample_used_in_res + how_many_samples_need_to_read < free);
	}
};

int main(int argc, char** argv)
{
	const float f_d = 7680000;

	bool isTx = false;
	float carrier_freq = 2.1276e9;
	float gain = 50;
	const float alpha = 0.22;

	int n;
	lms_info_str_t list[8]; //should be large enough to hold all detected devices
	if ((n = LMS_GetDeviceList(list)) < 0) error();
		
	lms_device_t *dev;
	if (LMS_Open(&dev, list[0], NULL)) error();
	if (LMS_Init(dev) != 0) error();
	//if (LMS_Calibrate(dev,isTx,0,1e6,0) != 0) error();
	if (LMS_SetAntenna(dev, isTx, 0, LMS_PATH_LNAW) != 0) error();
	if (LMS_EnableChannel(dev, isTx, 0, true) != 0) error();
	if (LMS_SetLOFrequency(dev, isTx, 0, carrier_freq) != 0) error();
	if (LMS_SetSampleRate(dev, f_d, 1) != 0) error();
	if (LMS_SetGaindB(dev, isTx, 0, gain) != 0) error(); //TX

	lms_stream_t sdr_stream; //stream structure
	sdr_stream.channel = 0; //channel number
	sdr_stream.fifoSize = 1024 * 1024 * 100; //fifo size in samples
	sdr_stream.throughputVsLatency = 1.0; //optimize for max throughput
	sdr_stream.isTx = isTx; //RX channel
	sdr_stream.dataFmt = lms_stream_t::LMS_FMT_I12; //12-bit integers

	ofstream out("dll.txt");
	wcdma_receiver receiver(1, f_d, 200, 100, 1.f, 0.5f,out);
	int sps = 2, sps_log2 = 1;
	int N = 1 << (sps_log2 + 16);
	std::vector<Ipp16sc> sig(N);
	if (LMS_SetupStream(dev, &sdr_stream) != 0) error();
	int receiver_needed_size;
	
	bool slow = true;
	while (1)
	{
		if (slow)
		{
			if (LMS_StartStream(&sdr_stream) != 0) error();
			for (int i = 0; i < 10; i++)
				LMS_RecvStream(&sdr_stream, sig.data(), N, NULL, 10000);
			LMS_StopStream(&sdr_stream);
			receiver_needed_size = receiver.general_work(sig,slow);
			if (LMS_StartStream(&sdr_stream) != 0) error();
			std::cout << "start demodulation" << std::endl;
			slow = false;
		}
		else
		{
			LMS_RecvStream(&sdr_stream, sig.data(), receiver_needed_size, NULL, 1000);
			receiver_needed_size = receiver.general_work(sig,slow);
		}
	}

	LMS_StopStream(&sdr_stream);
	LMS_DestroyStream(dev, &sdr_stream);
	LMS_Close(dev);

	_getch();
	return 0;	

}