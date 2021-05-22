# DSP-projects

1. IS-95 CDMA - Downlink Channel Model of IS-95 CDMA Standart written in LabVIEW 2017. Uses two NI USRP SDRs for real data 
transfer. Has a problem with linking to Shaping filter.dll.

2. FM Receiver - Digital Implementation of UHV analog radio with frequency modulation. 
Uses USRP SDR for listen to real radiostations.

3. filtration of markov chains - model for the nonlinear filtration algorithm  for binary Markov chains under gaussian noise.

4. WCDMA Receiver - model in LabVIEW 2017 for WCDMA Downlink Channels. It has a piece of real signal saved directly in test.vi, and
you can see correlation peaks for PSC & CPICH channels. Also there is a cpp real-time implementation of wcdma demodulator 
for scanning all possible channels. Demodulator needs Intel IPP library & LimeSuite API (now I can get signal only from 
LimeSDR).

5. DVB-T2 receiver - maybe implement in future.