# Fourier-and-Hilbert-transforms
Tools to analyze signals

There are 3 files at the moment.

- FFT.cpp is an implementation of the Fast Fourier transform (still missing the inverse transform implementation) with an example
in the main function which write out a file with the transformed vector of data.

- hilbert c++.cpp is an implementation of the Hilbert transform (using DFT and not FFT, still trying to implement FFT here and not DFT).
- modulation.cpp is an implementation of AM-modulation for signal including the SSB modulation and their respective demodulation. Notice
that in this last file I use functions defined in hilbert c++.cpp so it has to be included as header file .h if one want to use the functions
defined there.
