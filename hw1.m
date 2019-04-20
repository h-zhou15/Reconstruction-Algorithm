
clc,clear all;

readrawdata1;

%P(detector_Colsnum,detector_Rawsnum,beta_num)
%P(a,b,beta);
%%
N=256;
N_d=256;
beta_num=360;
pixelsize=2/N;  %2mm的圆
SOD=4;%光源到旋转中心的距离
SDD=8;%光源到探测器中心的距离
detectot_length=4;%探测器的长度
detector_channel_size=detectot_length/N_d;%探测器尺寸
fh_RL=zeros(1,2*N_d-1); %RL滤波函数

%--------------------------------
pro_beta=zeros(N_d,N_d);
weight_project_beta=zeros(N_d,N_d);
filtered_projection=zeros(N_d,N_d);
rec=zeros(N,N,N);
%%
%------------RL滤波核设计------%
for i=1:2*N_d-1
    fh_RL(i)=-1/(2*pi*pi*((i-N_d)*detector_channel_size).^2);
    if mod(i-N_d,2)==0
        fh_RL(i)=0;
    end
end
    fh_RL(N_d)=1/(8*detector_channel_size.^2);
    fh_RL=fh_RL.*detector_channel_size;
 %-------------------------------------------------
 %%
 tic
 for i=1:beta_num
     beta=(i-1)*pi/180;
     pro_beta=squeeze(p(:,:,i));
%      pro_beta=reshape(pro_beta,[N_d,N_d]);
     %========加权========%
     weight_project_beta=funcWeightProjectData(pro_beta,N_d,SDD,detector_channel_size);
     %=======卷积滤波====%
     filtered_projection=funcFilter(weight_project_beta,fh_RL,N_d);
     %=======反投影======%
     rec=rec+funcBackprojectionFdk(filtered_projection,pixelsize,SOD,beta,beta_num,N,detector_channel_size,SDD);
     t=toc;
     fprintf('theta: %d, timecost: %d s\n',i-1,t);
 end
    rec=rec*SDD/SOD;    %从虚拟探测器映射到实际探测器
	
    save('rec.mat','rec');
    toc
    imtool(squeeze(rec(:,:,N/2)));