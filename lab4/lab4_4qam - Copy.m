for snr_step = 1:1:10

    Nf = 1000; % in words
    % Na = 1000;(in bits)
    Na = 1000/2; % message word length (in symbel)
    T = 0.01; % T = 1ms
    eta = 64;
    fc = 800; %in Hz
    Eb = 1;
    Ts = T/eta;
    snr_var = snr_step; % in dB
    N0 = Eb * 10^(-snr_var/10);
    Nt = 100;
    phi_c = 2*pi*rand(1,1);


    v_train = 1 - 2*randi([0, 1], 1, Nt);
    totalErr = 0;
    q_sum = 0;
    for g = 1:Nf
        % data source
        a = randi([0 3], 1, Na);
        ah = zeros(1,Na);

        % symbol mapper = [+1 -1] BPSK
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
        v = [v_train v];
        Ns = eta*(Na+Nt);
        t = 0:T/eta:(Na+Nt)*T-T/eta;
        ht = (1/sqrt(T))*ones(1,eta);
        vt = conv(ht, upsample(v,eta));
        vt = vt(1:Ns);

        % modulator
        vct = real(vt * sqrt(2).*exp(1j*2*pi*fc*t));

        % channel
        rct = vct + sqrt(1/Ts * N0/2) * randn(1, length(vct));

        % demodulator
        rot = rct * sqrt(2).*exp(-1j*2*pi*fc*t);
        rot = rot * exp(1j*phi_c);

        % detector
        hrt = fliplr(ht);
        rt = conv(hrt,rot)*Ts;

        % decision device
        rn = rt(eta:eta:end);
        rn = rn/exp(1j*phi_c);
        rn_sync = rn(1:Nt);
        for u = 1:Nt
            q_sum = q_sum + rn_sync(u)/v_train(u);
        end
        q = q_sum/Nt;
        zn = rn(Nt+1:end)/q;
        for j = 1:Na
            if (real(zn(1,j)) >= 0) && (imag(zn(1,j)) > 0)
                ah(1,j) = 0;
            elseif (real(zn(1,j)) >= 0) && (imag(zn(1,j)) <= 0)
                ah(1,j) = 1;
            elseif (real(zn(1,j)) < 0) && (imag(zn(1,j)) <= 0)
                ah(1,j) = 3;
            elseif (real(zn(1,j)) <= 0) && (imag(zn(1,j)) >= 0)
                ah(1,j) = 2;
            end
        end


        % data sink
        totalErr = totalErr + biterr(a,ah);
    end

    % probability of error (theoretical)
    snr_linear = 10^(snr_var/10);
    Pb_thry = (1/2) * erfc(sqrt(snr_linear)); % qfunc(sqart(2)) = 1/2*erfc()

    % probability of error (calculated)
    Pb_calc = totalErr/(Na*2*Nf);

    Pb_thry_array(snr_var) = Pb_thry;
    Pb_calc_array(snr_var) = Pb_calc;
    % fprintf("error: %d\n",totalErr);

end

%plot
semilogy((1:snr_step),Pb_thry_array);
hold on;
semilogy((1:snr_step),Pb_calc_array);
hold off;
legend('Pb-thry','Pb-calc');
