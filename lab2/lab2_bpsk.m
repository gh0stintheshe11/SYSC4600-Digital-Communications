for snr_step = 1:1:10

    Nf = 1000; % in words
    Na = 1000; % (in bits)
    % Na = 1000/2; % message word length (in symbel)
    T = 0.01; % T = 1ms
    eta = 64;
    fc = 400; %in Hz
    % E[Vn] = 0.5 * (+1) + 0.5 * (-1) = 0
    % E[Vn^2] = 0.5 * (+1)^2 + 0.5 * (+1)^2 - E[Vn] = 1
    Eb = 1;
    Ts = T/eta;
    snr_var = snr_step; % in dB
    N0 = Eb * 10^(-snr_var/10);

    totalErrs = 0;
    for g = 1:Nf
        % data source
        a = randi([0 1], 1, Na);
        ah = zeros(1,Na);
        %fprintf('a : %s\n', sprintf('%d ', a));

        % symbol mapper = [+1 -1] BPSK
        % [A+Aj A-Aj -A+Aj -A-Aj] 4-QAM
        v = zeros(1,Na);
        for i = 1:Na
            if a(1,i) == 0
                v(1,i) = -1;
            elseif a(1,i) == 1
                v(1,i) = +1;
            end
        end

        % transmit filter
        Ns = eta*Na;
        t = 0:T/eta:Na*T-T/eta; %???
        ht = (1/sqrt(T))*ones(1,eta);
        vt = conv(ht, upsample(v,eta));
        vt = vt(1:Ns);

        % modulator
        vct = sqrt(2)*vt.*cos(2*pi*fc*t); %BPSK
        % vct = real(vt * sqrt(2).*exp(1j*2*pi*fc*t)); % 4-QAM

        % channel
        rct = vct + sqrt(1/Ts * N0/2) * randn(1, length(vct));

        % demodulator
        rot = sqrt(2)*rct.*cos(2*pi*fc*t); % BPSK
        % rot = rct * sqrt(2).*exp(1j*2*pi*fc*t); % 4-QAM
        rot = vt + sqrt(1/Ts * N0/2) * randn(1, length(vt));

        % detector
        hrt = fliplr(ht);
        rt = conv(hrt,rot)*Ts;

        % decision device
        vn = downsample(rt,eta);
        vn = vn(2:end);
        for j = 1:Na
            if vn(1,j) < 0
                ah(1,j) = 0;
            elseif vn(1,j) >= 0
                ah(1,j) = 1;
            end
        end

        % data sink
        %fprintf('a^: %s\n', sprintf('%d ',ah));
        %nErrs = biterr(a,ah);
        nErrs = sum(xor(a,ah));
        totalErrs = totalErrs + nErrs;
        %fprintf('Number of errors: %d\n', totalErrs);
    end
    % probability of error (theoretical)
    snr_linear = 10^(snr_var/10);
    Pb_thry = (1/2) * erfc(sqrt(snr_linear)); % qfunc(sqart(2)) = 1/2*erfc()
    %fprintf('Probability of errors (theoretical): %d\n', Pb_thry);

    % probability of error (calculated)
    Pb_calc = totalErrs/(Na*Nf);
    %fprintf('Probability of errors (calculated): %d\n', Pb_calc);

    Pb_thry_array(snr_var) = Pb_thry;
    Pb_calc_array(snr_var) = Pb_calc;


end

%plot
%hold on;

semilogy((1:snr_step),Pb_thry_array);
hold on
semilogy((1:snr_step),Pb_calc_array);
hold off;
legend('Pb-thry','Pb-calc');
