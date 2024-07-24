for snr_step = 1:1:10

    Nf = 1000; % in words
    % Na = 1000;(in bits)
    Na = 1000/2; % message word length (in symbel)
    T = 0.007; % T = 1ms
    eta = 64;
    L = 4;
    beta = 0.5;
    fc = 800; %in Hz
    Eb = 1;
    Ts = T/eta;
    snr_var = snr_step; % in dB
    No = Eb * 10^(-snr_var/10);

    totalErr = 0;
    PSD = 0;
    for g = 1:Nf
        % data source
        a = randi([0 3], 1, Na);
        ah = zeros(1,Na);

        % symbol mapper
        % [A+Aj A-Aj -A+Aj -A-Aj] 4-QAM
        v = zeros(1,Na);
        for i = 1:Na
            if a(1,i) == 0
                v(1,i) = 1+1j;
            elseif a(1,i) == 1
                v(1,i) = 1-1j;
            elseif a(1,i) == 2
                v(1,i) = -1+1j;
            elseif a(1,i) == 3
                v(1,i) = -1-1j;
            end
        end

        % transmit filter
        Ns = eta*Na;
        % t = 0:T/eta:Na*T-T/eta;
        hT = root_raised_cosine(beta, L, T, eta);
        vt = conv(hT, upsample(v,eta));
        vt = vt(1:L*eta+(Na-1)*eta);
        t = 0:T/eta:L*T+(Na-1)*T-T/eta;

        % modulator
        vct = real(vt * sqrt(2).*exp(1j*2*pi*fc*t));

        % channel
        rct = band_limited_channel(vct, Ts);
        rct = rct + sqrt(1/Ts*No/2)*randn(1, length(rct));

        % demodulator
        rot = rct * sqrt(2).*exp(-1j*2*pi*fc*t);

        % detector
        hrt = fliplr(hT);
        rt = conv(hrt,rot)*Ts;

        % decision device
        vn = rt(L*eta:eta:end);
        for j = 1:Na
            if (real(vn(1,j)) >= 0) && (imag(vn(1,j)) > 0)
                ah(1,j) = 0;
            elseif (real(vn(1,j)) >= 0) && (imag(vn(1,j)) <= 0)
                ah(1,j) = 1;
            elseif (real(vn(1,j)) < 0) && (imag(vn(1,j)) <= 0)
                ah(1,j) = 3;
            elseif (real(vn(1,j)) <= 0) && (imag(vn(1,j)) >= 0)
                ah(1,j) = 2;
            end
        end

        % data sink
        totalErr = totalErr + biterr(a,ah);
        %{
        Ns = length(vct); % Length of bandpass signal
        Vcf = fftshift(fft(vct)); % Calculate FFT
        PSD = Vcf.*conj(Vcf) * Ts / Ns .* sinc((-Ns/2:Ns/2-1)/Ns).^2 + PSD;
        %}
    end
    % probability of error (theoretical)
    snr_linear = 10^(snr_var/10);
    Pb_thry = (1/2) * erfc(sqrt(snr_linear)); % qfunc(sqart(2)) = 1/2*erfc()

    % probability of error (calculated)
    Pb_calc = totalErr/(Na*2*Nf);

    Pb_thry_array(snr_var) = Pb_thry;
    Pb_calc_array(snr_var) = Pb_calc;

    %fprintf('error: %s\n', sprintf('%d ', totalErr));


end
%{
f = (-Ns/2:Ns/2-1) / (Ns*Ts); % Frequencies
nexttile;
plot(f, 10*log10(PSD/Nf));
axis([650 950 -80 10]);
%}

nexttile
semilogy((1:snr_step),Pb_thry_array);
hold on;
semilogy((1:snr_step),Pb_calc_array);
hold off;
legend('Pb-thry','Pb-calc');
