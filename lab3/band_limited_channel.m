function rct = band_limited_channel(vct, Ts)
%BAND_LIMITED_CHANNEL Simulates a bandpass channel with limited bandwidth.
%
%   rct = band_limited_channel(vct, Ts) simulates the transmission of
%   a bandpass signal, given by vct, sampled at a rate of Ts samples
%   per second, through a channel with a bandwidth of 150 Hz centered
%   around 800 Hz.
    fc = 800;                   % Center frequency
    BW = 250;                   % Channel bandwidth
    Ns = length(vct);           % Length of transmitted signal (samples)

    f = (-Ns/2:Ns/2-1) / (Ns*Ts);

    % Channel frequency response
    H = zeros(1, Ns);
    H(abs(abs(f)-fc) < BW/2) = 1;

    Vcf = fftshift(fft(vct));
    Rcf = H .* Vcf;
    rct = ifft(ifftshift(Rcf));
end
