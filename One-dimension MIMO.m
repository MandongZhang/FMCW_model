%{
Data:2022/12/2
Function: Angle Measurement
Imprefection: 
    1. The algorithm Can not get the direction of velocity. And it just get
    the absolute speed.
    2. It can not display angle from -90 degrees to 90 degress.
    3. It is not sampling sufficiently.
%}

clc
close all
clear all
% parameter 
fc = 17e+7;
c = 3e+8;
B = 3e+7;
Ns = 2^10; %need to set
Nc = 256;
Tchirp = 4e-4;
t = (0:Ns*Nc-1)*Tchirp/Ns;
S = B/Tchirp; 

vres = c/(2*Nc*Tchirp*fc);
rangeRes = c/(2*B);
Fs = Ns/Tchirp;
R_max = 200;
f_IFmax = 2*S*R_max/c;

%Transmitter
lambda = c/fc;
N_Rx = 128;
Tx = cos(2*pi*fc.*(t)+pi*S.*mod(t, Tchirp).^2);
position = [0,0];
velocity = [0,0];
[R,Theta] = ObjectGenerate(t, position, velocity); % R and Theta are  0D values
tau = 2*R/c;
D_Rx = lambda/2;    %Distance of two adjacent receivers
D_phase = 2*pi*D_Rx*sind(Theta(1))/lambda; % The difference of phase between two adjacent receivers
angleRes =  asind(lambda/(N_Rx*D_Rx));

%Receiver
Rx = zeros( 1,Ns*Nc, N_Rx);
Mix = zeros(1, Ns*Nc, N_Rx);

for i=1:N_Rx
     Rx(1,:,i) = cos(2*pi*fc*(t-tau)+pi*S.*mod((t-tau), Tchirp).^2+(i-1)*D_phase); 
     Mix(1,:,i) = Rx(1,:,i).*Tx;
     %Mix(1,:,i) = lowpass(Mix(1,:,i),1.1*f_IFmax, Fs);
end

Mix = reshape(Mix,[Ns,Nc, N_Rx]); 
Data_VandR = Mix(:,:,1);

Data_RandA = zeros(Ns, N_Rx);
for j = 1:N_Rx
   Data_RandA(:,j) = Mix(:,1,j); 
end

VandR_fft = abs(fft2(Data_VandR));
VandR_fft = fftshift(VandR_fft);
VandR_fft = VandR_fft(Ns/2+1:Ns,:);
RandA_fft = abs(fft2(Data_RandA));
RandA_fft = fftshift(RandA_fft);
RandA_fft = RandA_fft(Ns/2+1:Ns,:);


figure(1)
n_range = linspace(0,Ns/2-1, Ns/2); % can not display -90 degrees.
range_axis = rangeRes*n_range; 
n_v = linspace(-Nc/2, Nc/2-1, Nc);
v_axis = vres*n_v;
mesh(v_axis,range_axis,VandR_fft);
title('range and velocity');
xlabel('velocity');
ylabel('range');


figure(2)
n_angle = linspace(N_Rx/2, -N_Rx/2+1, N_Rx); 
angle_axis = asind(n_angle*lambda/(N_Rx*D_Rx));
mesh(angle_axis,range_axis,RandA_fft);
title('range and angle');
xlabel('angle');
ylabel('range');






