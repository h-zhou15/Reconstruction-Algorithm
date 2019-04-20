function [Img] = fbp(N,projection,filtertype,detector_size,pixelsize)
    M=size(projection,1);
    N_d=M;
    TT=linspace(-N_d*detector_size/2,N_d*detector_size/2,N_d+1);    %探测器边界坐标
    for i=1:N_d
        DD(i)=(TT(i)+TT(i+1))/2;    %探测器通道的中点坐标
    end

    Img=zeros(N);
    Detectorsize=detector_size;
    [ImgRow,ImgCol]=size(Img);
    width=2^nextpow2(size(projection,1))*16;
    pro_fft=fft(projection,width);
    
    R_L =(1/Detectorsize)*[0:(width/2), width/2-1:-1:1]'/width; %1/2tau,tau为采样间距
    filter_Hamming=0.54-0.46*cos([-width/2+1:0,1:width/2]'.*(2*pi/(width)));
    pro_filter=zeros(width,size(projection,2));
    if filtertype=='None'
        for i=1:1:size(projection,2)
            %pro_filter(:,i)=pro_fft(:,i).*R_L;
            pro_filter(:,i)=pro_fft(:,i).*R_L.*filter_Hamming;
        end
    end
    pro_ifft=real(ifft(pro_filter));
    for i=1:1:size(projection,2)
        for c=1:ImgCol
            for r=1:ImgRow
                x=(c-(ImgCol+1)/2).*pixelsize;
                y=((ImgRow+1)/2-r).*pixelsize;
                s=x*cos((i-1)*pi/180)+y*sin((i-1)*pi/180);
                s=s/Detectorsize+N_d/2+1;
                n=floor(s);
                if n>0&&n<=N_d
                  q=pro_ifft(n,i);
                  Img(r,c)=Img(r,c)+q;
                end
            end
        end
    end
    Img=Img*pi/size(projection,2);
%     
%     %cbp
%     %%=================R_L================%%
%     N_d=size(projection,1);
%     tau=Detectorsize;
%     fh_RL=zeros(2*N_d-1,1);
%     for k=1:2*N_d-1
%         fh_RL=-1/(((k-N_d)*pi*tau).^2);
%         if mod(k-N_d,2)==0
%             fh_RL(k)=0;
%         end
%     end
%     fh_RL(N_d)=1/(4*tau*tau);
%     fh_RL=fh_RL*tau;
%     K=size(projection,2);
%     delta=1;
%     for i=1:delta:K
%        P=projection(:,i);
%        k=N_d:2*N_d-1;
%        q=conv(P,fh_RL);
%        q=q(k);
%        cm=(N_d/2)*(1-cos((i-1)*pi/180)-sin((i-1)*pi/180));
%        for c=1:N
%            for r=1:N
%                x=(c-(ImgCol+1)/2);
%                y=((ImgRow+1)/2-r);
%                s=x*cos((i-1)*pi/180)+y*sin((i-1)*pi/180);
%                s=s/Detectorsize+N_d/2+1;
%                %n=floor(s); 
%                n=round(s);
%                if n>0&&n<=N_d
%                  Img(r,c)=Img(r,c)+q(n);
%                end
%            end
%        end
%     end
%     Img=Img*pi/size(projection,2);
    
end

