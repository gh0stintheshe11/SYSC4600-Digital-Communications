for snr_step = 1:1:14

    Nf = 1000; % in words
    % Na = 1000;(in bits)
    Na = 1000/4; % message word length (in symbel)
    T = 0.01; % T = 1ms
    eta = 64;
    fc = 400; %in Hz
    Eb = 1;
    Ts = T/eta;
    snr_var = snr_step; % in dB
    N0 = Eb * 10^(-snr_var/10);

    totalErr = 0;
    for g = 1:Nf
        % data source
        a = randi([0 15], 1, Na);
        ah = zeros(1,Na);

        % symbol mapper
        % [-3+3j -1+3j | +1+3j +3+3j]
        % [-3+1j -1+1j | +1+1j +3+1j]
        % ---------------------------
        % [-3-1j -1-1j | +1-1j +3-1j]
        % [-3-3j -1-3j | +1-3j +3-3j] 16-QAM
        v = zeros(1,Na);
        for i = 1:Na
            % 0~3
            if a(1,i) == 0
                v(1,i) = -3+3j;
            elseif a(1,i) == 1
                v(1,i) = -3+1j;
            elseif a(1,i) == 2
                v(1,i) = -3-1j;
            elseif a(1,i) == 3
                v(1,i) = -3-3j;
                % 4~7
            elseif a(1,i) == 4
                v(1,i) = -1+3j;
            elseif a(1,i) == 5
                v(1,i) = -1+1j;
            elseif a(1,i) == 6
                v(1,i) = -1-1j;
            elseif a(1,i) == 7
                v(1,i) = -1-3j;
                % 8~11
            elseif a(1,i) == 8
                v(1,i) = +1+3j;
            elseif a(1,i) == 9
                v(1,i) = +1+1j;
            elseif a(1,i) == 10
                v(1,i) = +1-1j;
            elseif a(1,i) == 11
                v(1,i) = +1-3j;
                % 12~15
            elseif a(1,i) == 12
                v(1,i) = +3+3j;
            elseif a(1,i) == 13
                v(1,i) = +3+1j;
            elseif a(1,i) == 14
                v(1,i) = +3-1j;
            elseif a(1,i) == 15
                v(1,i) = +3-3j;
            end
        end

        % transmit filter
        Ns = eta*Na;
        t = 0:T/eta:Na*T-T/eta; %???
        ht = (1/sqrt(T))*ones(1,eta);
        vt = conv(ht, upsample(v,eta));
        vt = vt(1:Ns);

        % modulator
        vct = real(vt * sqrt(2).*exp(1j*2*pi*fc*t)); % QAM

        % channel
        rct = vct + sqrt(1/Ts * N0/2) * randn(1, length(vct));

        % demodulator
        rot = rct * sqrt(2).*exp(-1j*2*pi*fc*t); % QAM

        % detector
        hrt = fliplr(ht);
        rt = conv(hrt,rot)*Ts;

        % decision device
        vn = rt(eta:eta:end);
        for j = 1:Na
            % 0
            if real(vn(1,j)) <-2 ...
                    && imag(vn(1,j)) >2
                ah(1,j) = 0;
                % 1
            elseif real(vn(1,j)) <-2 ...
                    && imag(vn(1,j)) <2 ...
                    && imag(vn(1,j)) >0
                ah(1,j) = 1;
                % 2
            elseif real(vn(1,j)) <-2 ...
                    && imag(vn(1,j)) <0 ...
                    && imag(vn(1,j)) >-2
                ah(1,j) = 2;
                % 3
            elseif real(vn(1,j)) <-2 ...
                    && imag(vn(1,j)) <-2
                ah(1,j) = 3;
                % 4
            elseif real(vn(1,j)) >-2 ...
                    && real(vn(1,j)) <0 ...
                    && imag(vn(1,j)) >2
                ah(1,j) = 4;
                % 5
            elseif real(vn(1,j)) >-2 ...
                    && real(vn(1,j)) <0 ...
                    && imag(vn(1,j)) <2 ...
                    && imag(vn(1,j)) >0
                ah(1,j) = 5;
                % 6
            elseif real(vn(1,j)) >-2 ...
                    && real(vn(1,j)) <0 ...
                    && imag(vn(1,j)) <0 ...
                    && imag(vn(1,j)) >-2
                ah(1,j) = 6;
                % 7
            elseif real(vn(1,j)) >-2 ...
                    && real(vn(1,j)) <0 ...
                    && imag(vn(1,j)) <-2
                ah(1,j) = 7;
                % 8
            elseif real(vn(1,j)) >0 ...
                    && real(vn(1,j)) <2 ...
                    && imag(vn(1,j)) >2
                ah(1,j) = 8;
                % 9
            elseif real(vn(1,j)) >0 ...
                    && real(vn(1,j)) <2 ...
                    && imag(vn(1,j)) <2 ...
                    && imag(vn(1,j)) >0
                ah(1,j) = 9;
                % 10
            elseif real(vn(1,j)) >0 ...
                    && real(vn(1,j)) <2 ...
                    && imag(vn(1,j)) <0 ...
                    && imag(vn(1,j)) >-2
                ah(1,j) = 10;
                % 11
            elseif real(vn(1,j)) >0 ...
                    && real(vn(1,j)) <2 ...
                    && imag(vn(1,j)) <-2
                ah(1,j) = 11;
                % 12
            elseif real(vn(1,j)) >2 ...
                    && imag(vn(1,j)) >2
                ah(1,j) = 12;
                % 13
            elseif real(vn(1,j)) >2 ...
                    && imag(vn(1,j)) <2 ...
                    && imag(vn(1,j)) >0
                ah(1,j) = 13;
                % 14
            elseif real(vn(1,j)) >2 ...
                    && imag(vn(1,j)) <0 ...
                    && imag(vn(1,j)) >-2
                ah(1,j) = 14;
                % 15
            elseif real(vn(1,j)) >2 ...
                    && imag(vn(1,j)) <-2
                ah(1,j) = 15;
            end
        end

        % data sink
        totalErr = totalErr + biterr(a,ah);
    end
    % probability of error (theoretical)
    snr_linear = 10^(snr_var/10);
    Pb_thry = (1/2) * erfc(sqrt(snr_linear)); % qfunc(sqart(2)) = 1/2*erfc()

    % probability of error (calculated)
    Pb_calc = totalErr/(Na*4*Nf);

    Pb_thry_array(snr_var) = Pb_thry;
    Pb_calc_array(snr_var) = Pb_calc;

end

%plot
semilogy((1:snr_step),Pb_thry_array);
hold on;
semilogy((1:snr_step),Pb_calc_array);
hold off;
legend('Pb-thry','Pb-calc');

