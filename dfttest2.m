clear
delete(instrfindall);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize communication with Arduino
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arduino = serial('COM5'); %fill in the blank with the serial port used
arduino.Baudrate=9600;
fopen(arduino);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize sample and counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleSize = 10500;
a=1;
dataECG = zeros(sampleSize,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Collection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (a <= sampleSize)
disp(a/sampleSize*100);
ecgRaw = fgets(arduino);
dataECG(a,1) = str2double(ecgRaw)/1023*5;
a = a + 1;
end
sfz = 180;
ecg = dataECG;
t = 1/sfz * (1:size(ecg,1));
secg = ecg(1000:4000,1); %stores the 1000th sample to 4000th sample from "ecg" in "secg"
fz = sfz/length(secg); %calculation of frequency "fz"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment 2
% Plot Fourier Transform of acquired data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fft_secg = zeros(length(secg),1); %initializes 2D array "fft_secg" with same length as "secg"
disp('Calculating Fourier Transform');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Insert Fourier Transform Code below%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(secg);
sum = 0;
for k=1:N
    for n=1:N
        add = secg(n)*exp(-(1j*2*pi*(k-1)*(n-1))/N);
        sum = sum + add;
    end
    fft_secg(k) = sum;
    sum = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Insert Fourier Transform Code above%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Done');

fsaxis = (0:floor(length(secg)/2))*fz; %mirroring of frequency domain
fsaxis = [fsaxis -fliplr(fsaxis(2:end))]; %ensure no spacing between “-“ and “fliplr” function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical Display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
%subplot(121);
plot(t, ecg(:,1));
xlabel('time (sec)');
ylabel('x(t) magnitude');
%subplot(122);
figure;
plot(fsaxis(2:end), abs(fft_secg(2:end)));
xlabel('frequency (Hz)');
ylabel('X(jw) magnitude');
drawnow;

