% Question 1: xor(a,ah) is doing a bitwise xor between a and ah, so if ah
% is same as a then 0, not same yells 1, sum all 1s = sum all differet
% sents and recives = sum all errors

Nf = 10;

Na = 128; % message word length (in bits)
T = 0.01; % T = 1ms
eta = 64;
fc = 400;

    % data source
    a = randi([0 1], 1, Na);
    % a = [1 0 1 0 1 0 1 1 0 1];
    ah = zeros(1,Na);
    fprintf('a : %s\n', sprintf('%d ', a));

    % symbol mapper
    v = zeros(1,Na);
    for i = 1:Na
        if a(1,i) == 0
            v(1,i) = 1;
        elseif a(1,i) == 1
            v(1,i) = -1;
        end
    end

    % transmit filter
    Ns = eta*Na;
    t = 0:T/eta:Na*T-T/eta; %???
    ht = (1/sqrt(T))*ones(1,eta);
    vt = conv(ht, upsample(v,eta));
    vt = vt(1:Ns);
    tiledlayout(3,2);
    nexttile;
    plot(vt);
    title('v(t)');
    ylim([-1.1, 1.1]/sqrt(T));

    %modulator
    vct = sqrt(2)*vt.*cos(2*pi*fc*t);
    nexttile;
    plot(vct);
    title('vc(t)');

    %ideal digital channel
    rct = vct;

    % demodulator
    rot = sqrt(2)*rct.*cos(2*pi*fc*t);
    nexttile;
    plot(rot);
    title('ro(t)');

    % detector
    Ts = T/eta;
    hrt = fliplr(ht);
    rt = conv(hrt,rot)*Ts;
    nexttile;
    plot(rt);
    title('r(t)');
    ylim([-1.1, 1.1]);

    % decision device
    vn = downsample(rt,eta);
    vn = vn(2:end);
    for j = 1:Na
        if vn(1,j) < 0
            ah(1,j) = 1;
        elseif vn(1,j) >= 0
            ah(1,j) = 0;
        end
    end


% data sink
fprintf('a^: %s\n', sprintf('%d ',ah));
nErrs = sum(xor(a,ah));
fprintf('Number of errors: %d\n', nErrs);

% step 5
Vf = fftshift(fft(vt));
freq = (-Ns/2:Ns/2-1)/(Ns*Ts);
PSD = abs(Vf).^2 * Ts / Ns .* sinc((-Ns/2:Ns/2-1)/Ns).^2;
vt_db = 10*log10(PSD);
PSD_thry = sinc(freq*T).^2;
PSD_thry_db = 10*log10(PSD_thry);
nexttile;
hold on;
plot(freq,vt_db);
plot(freq,PSD_thry_db);
title('PSD');
legend('PSD','PSD_thry');
xlim([-750 750]);
ylim([-50 10]);


