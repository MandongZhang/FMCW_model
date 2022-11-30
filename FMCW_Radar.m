%{
Data:2022/11/28
Author:Mandong Zhang
Name:FMCW Radar System Model
Function:Generateing several moving objects,detecting the range and
velocity
Imperfection:
    1.If the range is too close, aliasing will happen because the
sampling after Mixer.
    2.Direct compensation of range is not suitable for multiple objects.
%}

clc
clear

fc = 17e+7;
c = 3e8;
Tchirp = 4e-4;
B = 3e7;
Nc = 2^9;
Ns = ceil(6*(fc+B)*Tchirp);

%Ŀ������
%R_init = [100;100]; 
%v = [30;50];
R_init = 10;
v = 15;
n = 0:Ns*Nc-1;
t = n.*(Tchirp/Ns);

%��֤����
S  = B/Tchirp; 
rangeRes = c/(2*B);
vres = c/(2*Nc*Tchirp*fc);
Fs = Ns/Tchirp;
f_IF = 2*S*R_init./c;
R_exact = R_init+v*Tchirp*Nc; %ʵ�������λ��
R_max = 200;
f_IFmax = 2*S*R_max/c;

%�����ź�����
Tx = cos(2*pi*fc.*(t)+pi*S.*mod(t, Tchirp).^2);

%�ز�����
tau = 2*(R_init+v*t)./(c+v);
Rx = cos(2*pi*fc.*(t-tau) + S*pi.*mod((t-tau), Tchirp).^2);
%Rx = sum(Rx);

%��Ƶ �˲� 
Mix = Tx.*Rx;
Mix = lowpass(Mix,1.5*f_IFmax, Fs);
Mix = reshape(Mix, Ns, Nc);

%��ȡ���������fft
N_abstract = 2^10;
Abstract_Interval = floor(Ns/N_abstract);
Abstract_Series = 1:Abstract_Interval:Abstract_Interval*N_abstract;
F_abstract = length(Abstract_Series)/Tchirp;
Mix = Mix(Abstract_Series, :);

%��άfft�����ٶȾ���
sig_fft = abs(fft2(Mix));
sig_fft = fftshift(sig_fft);
sig_fft = sig_fft(N_abstract/2:N_abstract-1,:);

n_range = linspace(0,N_abstract/2-1, N_abstract/2);
range_axis = rangeRes*n_range; 
n_v = linspace(-Nc/2, Nc/2-1, Nc);
v_axis = vres*n_v;
mesh(v_axis,range_axis,sig_fft);

[row, col] = find(sig_fft==MaxValue);
R_solve = range_axis(row);
v_solve = v_axis(col);
R_comp = R_solve - fc*v_solve/S;
fprintf('R = %d, v = %d', R_comp,v_solve);




