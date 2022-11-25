
%% 公众号 ：调皮连续波
%% 时间：2022年10月23日
clc;
clear all;
close all;

Frame =1; %帧数设置；

for frame =1:Frame %帧数设置

%% 雷达参数设置
parameter  = generateParameter();
parameter.frame = frame;
%% 雷达回波信号建模
rawData    = generateSignal(parameter);
firstChirp = rawData(1,:,1);

% figure(1);
% plot(real(firstChirp));
% hold on;
% plot(imag(firstChirp));
% xlabel('采样点数'); ylabel('幅值');title('原始数据实虚部');%第1个chirp。

%% IQ通道校正
firstChirp_I = imag(firstChirp);
firstChirp_Q = real(firstChirp);
%设定幅度误差因子
alpha = 2;
%设定相位误差因子
phi = 5;
firstChirp_Q_1 = (alpha+1)*firstChirp_Q*exp(phi*2*pi/360);
%% IQ通道不平衡 效果
firstChirp_IQ_1 = complex(firstChirp_Q_1,firstChirp_I);
%IQ通道校正算法
%求均值
firstChirp_I = imag(firstChirp_IQ_1);
firstChirp_Q = real(firstChirp_IQ_1);
I_before_correction =firstChirp_I-mean(firstChirp_I);
Q_before_correction=firstChirp_Q -mean(firstChirp_Q);
%估计参数
alphi = sqrt(mean(Q_before_correction.*Q_before_correction)/mean(I_before_correction.*I_before_correction))-1;
phi = -asin(mean(I_before_correction.*Q_before_correction)/sqrt(mean(I_before_correction.*I_before_correction)*mean(Q_before_correction.*Q_before_correction)));
%P矩阵求解
P = [1,0;tan(phi),1/((1+alphi)*cos(phi))];
%计算IQ
IQ = P*[I_before_correction;Q_before_correction];

%重组信号
I =IQ(1,:);
Q =IQ(2,:);
signal_IQ =Q+I*1j;

%图形绘制
% figure(2)
% subplot(2,1,1);
% plot((abs(fft(firstChirp_IQ_1))));
% xlabel('距离(m)'); ylabel('幅值');title('距离维FFT(校正前)');
% subplot(2,1,2);
% plot((abs(fft(signal_IQ))));
% xlabel('距离(m)'); ylabel('幅值');title('距离维FFT(校正后)');
%%

%% 雷达信号处理
rangeRes     = parameter.c / (2 * parameter.BandwidthValid); %距离分辨率 有效带宽
rangeIndex   = (0:parameter.rangeBin-1) * rangeRes;
speedRes     = parameter.lambda / (2 * parameter.dopplerBin * parameter.Tr);
dopplerIndex = (-parameter.dopplerBin/2:1:parameter.dopplerBin/2 - 1) * speedRes;
angleRes     = parameter.lambda / (parameter.virtualAntenna * parameter.dx) * 180 / pi;
angleIndex   = (-parameter.virtualAntenna/2:1:parameter.virtualAntenna/2 - 1) * angleRes;

%%1D FFT
fft1dData    = fft(firstChirp);
% figure(3);
% plot(rangeIndex,db(abs(fft1dData)./max(abs(fft1dData))));
% xlabel('距离(m)'); ylabel('幅值(dB)');title('距离维FFT');
% 

%% 2D FFT
%% 距离-多普勒谱
channelNum    = size(rawData,1);
rangebinNum   = size(rawData,2);
dopplerbinNum = size(rawData,3);
fft2dDataPower= zeros(size(rawData));
fft2dDataDB   = zeros(size(rawData));
fftRADataPower= zeros(size(rawData));
for chanId = 1:1:channelNum
    fft2dDataPower(chanId,:,:) = RDfftMatrix(rawData(chanId,:,:));
end

% figure(4);
% mesh(dopplerIndex',rangeIndex,db(abs(squeeze(fft2dDataPower(chanId,:,:)))));
% view(2);
% xlabel('速度(m/s)'); ylabel('距离(m)'); zlabel('幅值');
% title('距离-多普勒谱');
% mesh(abs(squeeze(fft2dDataPower(chanId,:,:))));

%% 距离-角度谱
for dopplerId = 1:1:dopplerbinNum
    fftRADataPower(:,:,dopplerId) = RAfftMatrix(rawData(:,:,dopplerId));
end
% figure(5);
% mesh(rangeIndex,angleIndex,(abs(squeeze(fftRADataPower(:,:,dopplerId)))));
% view(2);
% xlabel('距离(m)'); ylabel('角度'); zlabel('幅值');
% title('距离-角度谱');

% RAMap = fftRADataPower(:,:,1);
% for angle =1:12
%     for range =1:rangebinNum
%      X(angle,range) = RaMap(angle,range)*cos( RAMap(angle,range));
%      Y(angle,range)= RaMap(angle,range)*sin( RAMap(angle,range));
%     end
% end
% figure();
% mesh(rangeIndex,angleIndex,abs(Y));


%% 多通道非相干积累
accumulateRD = chan_Accumulate((fft2dDataPower));
% mesh(dopplerIndex',rangeIndex,(accumulateRD));
% figure(1);
% mesh(dopplerIndex',rangeIndex,db(accumulateRD));
% view(2);
% xlabel('速度(m/s)'); ylabel('距离(m)'); zlabel('幅值');
% title(['通道积累 第',num2str(frame),'帧']);
% pause(0.01);

%% CFAR检测
cfarParameter = generateCfarParameter(); %生成cfar数据
[pointList,cfarRD] = cfar(cfarParameter,db(accumulateRD));
% figure;
% mesh(dopplerIndex',rangeIndex,cfarRD);
% xlabel('速度(m/s)'); ylabel('距离(m)'); zlabel('幅值');
% title('cfar');

%% peakSearch
[RD_pearkSearch,peakSearchList] = peakSearch(cfarRD,pointList);
% figure;
% mesh(dopplerIndex',rangeIndex,RD_pearkSearch);
% xlabel('速度(m/s)'); ylabel('距离(m)'); zlabel('幅值');
% title('峰值搜索');
detectPointNum = size(peakSearchList,2);

%% DOA估计
for targetIdx = 1:detectPointNum
    rangeBin = peakSearchList(1,targetIdx);
    speedBin = peakSearchList(2,targetIdx);
    range = (rangeBin - 1) * rangeRes;
    speed = (speedBin - parameter.dopplerBin/2 - 1) * speedRes;

    ant = squeeze(fft2dDataPower(:,rangeBin,speedBin));
    compCoff = generateCompCoff(parameter,speedBin); %多普勒补偿系数
  
    if 1 %是否进行多普勒补偿 
        ant = ant .* compCoff;
        [angle,doa_abs] = doa(parameter,ant);
    end

    figure(6);
    angleIndex =((-512:1:512-1)/512) * 180 / pi;
    hold on;
    plot(angleIndex,doa_abs);
    title('多普勒补偿前后测角结果');
    xlabel('角度');ylabel('幅值');
    fprintf('目标%d的距离为%f,速度为%f,角度为%f\n',targetIdx,range,speed,angle);
end

end
