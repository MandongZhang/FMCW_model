%{
Data:2022/12/2
Function: Angle Measurement
Imprefection: 
    1. Can not distingush postive and negetive angle.
    2. X axis can not display the range of angle above 80 degree.
%}

fc = 17e+7;
c = 3e+8;
Ns = 2^10;
Nc = 8;
t = 0:Ns*Nc-1;

%Transmitter
lambda = c/fc;
Tx = cos(2*pi*fc*t);

%Object parameter
Theta = 90;
R = 400;
tau = 2*R/c;
N_Rx = 128;

%Receiver
Rx = zeros(N_Rx, Ns*Nc);
D_Rx = lambda/2;    %Distance of two adjacent receivers
D_phase = 2*pi*D_Rx*sin(2*pi*Theta/360)/lambda; % The difference of phase between two adjacent receivers


for i=1:N_Rx
    Rx(i,:) = cos(2*pi*fc*(t-tau)+i*D_phase);
end

Rx_fft = abs(fft(Rx));
angle_axis = asin((0:N_Rx/2-1)*lambda/(N_Rx*D_Rx))*180/pi;
figure(1)
mesh(1:Ns*Nc,angle_axis,Rx_fft(1:N_Rx/2,:))
title('Angle Measurement theta = 80бу');
ylabel('angle positive');
figure(2)
mesh(Rx_fft);