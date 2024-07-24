%
% Data Source
%
Na = 128; % Message word length (in bits)
T=0.01; %transfer ms to second
eta=64;
%a = [1,0,1,0,1,0,1,1,0,1]; %for step 1 to 4
a = randi([0 1], 1, Na); %for step 5
fc=400;
Ts=T/eta;
Ns=eta*Na;
disp(sprintf('Transmitted Message: %s', sprintf('%d ', a)));
%
%Symbol mapper
%
for i=1 : Na
    if a(1,i)==1
        v(1,i)=-1;
    else v(1,i)=1;
    end
end
%
%transimit filter
%
t=0:T/eta:Na*T-T/eta; %set up time period for hT
hT=1/sqrt(T)*ones(1,eta); %Sampled Normalized pulse
vt=conv(upsample(v,eta),hT);
vt=vt(1:(eta*Na));
%
%Modulator
%
vc=sqrt(2)*vt.*cos(2*pi*fc*t);
%
% Ideal Digital Channel
%
rc=vc;
rn = v; % for step2 only
%
%Demodulator
%
r0 = sqrt(2)*rc.*cos(2*pi*fc*t);
%
%Detector
%
hrT=fliplr(hT);
rt=conv(hrT,r0)*(T/eta);
rnt=rt(eta:eta:eta*Na);
%Decision Device
for i=1 : Na
    if rnt(1,i) >= 0 
        ah(1,i)=0; 
    else ah(1,i)=1; 
    end 
end
%
% Data Sink
%

%plot(t,vt); %for question 2 only
%ylim([-1.1/sqrt(T) 1.1/sqrt(T)]); %for question 2 only
Vf = fftshift(fft(vt));
freq=(-Ns/2:Ns/2-1)/(Ns*Ts);
PSD = abs(Vf).^2 * Ts / Ns .* sinc((-Ns/2:Ns/2-1)/Ns).^2;
dBvt=10*log10(PSD); %convert experimental PSD in dB unit
theoreyPSD= sinc(freq*T).^2;
dBtheorey=10*log10(theoreyPSD);
hold on
plot(freq,dBvt);
plot(freq, dBtheorey);
ylim([-50 10]);% for question 7
xlim([-750 750]); % for question 7
legend('Experimental','Theoretical');
hold off
disp(sprintf('Received Message: %s', sprintf('%d ', ah)));
nErrs = sum(xor(a, ah));
disp(sprintf('Number of errors: %d', nErrs));
