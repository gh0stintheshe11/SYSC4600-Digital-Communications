%
% Data Source
%
Na = 128; % Message word length (in bits)
Nf = 1000; %question 9
T=0.01; %transfer ms to second
eta=64;
fc=400;
Ts=T/eta;
Ns=eta*Na;
a=zeros(Nf,Na);
v=zeros(Nf,Na);
tot=zeros(1,Ns);
Vf=zeros(Nf,Ns);
%a = [1,0,1,0,1,0,1,1,0,1]; %for step 1 to 4
for i=1:Nf
    
    a(i,1:Na) = randi([0 1], 1, Na); %for step 5

    %disp(sprintf('Transmitted Message: %s', sprintf('%d ', a)));
%
%Symbol mapper
%
    for x=1 : Na
        if a(i,x)==1
            v(i,x)=-1;
        else 
            v(i,x)=1;
        end
    end
%
%transimit filter
%

    t=0:T/eta:Na*T-T/eta; %set up time period for hT

    hT=1/sqrt(T)*ones(1,eta); %Sampled Normalized pulse
    vt(i,:)=conv(upsample(v(i,1:Na),eta),hT);
    vt1(i,1:(eta*Na))=vt(i,1:(eta*Na));
%
%Modulator
%
    vc=sqrt(2)*vt1(i,:).*cos(2*pi*fc*t);
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
    for y=1 : Na
        if rnt(1,y) >= 0 
            ah(i,y)=0; 
        else ah(i,y)=1; 
        end 
    end
%
% Data Sink
%

%plot(t,vt); %for question 2 only
%ylim([-1.1/sqrt(T) 1.1/sqrt(T)]); %for question 2 only


Vf(i,1:Ns) = fftshift(fft(vt1(i,:)));
freq=(-Ns/2:Ns/2-1)/(Ns*Ts);
PSD(i,1:Ns) = abs(Vf(i,1:8192)).^2 * Ts / Ns .* sinc((-Ns/2:Ns/2-1)/Ns).^2;

%dBvt=10*log10(PSD(i)); %convert experimental PSD into dB unit
tot=(tot+PSD(i,:));
end
avg=tot./Nf;
theoreyPSD= (sinc(freq*T)).^2; %the equation of theroyPSD
dBtheorey=10*log10(theoreyPSD);%convert theroy PSD into dB unit
hold on
plot(freq,10*log10(avg));%question 7,8,9 plot
plot(freq, dBtheorey);%question 7,8,9 plot
ylim([-50 10]);% for question 7,8,9
xlim([-750 750]); % for question 7,8,9
legend('Experimental','Theoretical');
hold off

%disp(sprintf('Received Message: %s', sprintf('%d ', ah)));
nErrs = sum(xor(a, ah));
%disp(sprintf('Number of errors: %d', nErrs));
