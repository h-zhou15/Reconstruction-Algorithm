clc,clear;

M=512;%探测器的个数
A=200;C=150;
B=sqrt(A^2-C^2);

%首先画出椭圆的图像
Img=zeros(512);

for y=1:M
    for x = 1:M
        dis=sqrt((x-(M/2-C))^2+(y-M/2)^2)+sqrt((x-(M/2+C))^2+(y-M/2)^2);    %竖直方向为x轴，img
        
        if dis<=2*A
            Img(x,y)=10;
        else
            Img(x,y)=0;
        end
    end
end


%投影值
%projection=zeros(512,180);

for i=1:180
    theta=i*pi/180;
    r2=A^2*(cos(theta))^2+B^2*(sin(theta))^2;
    for s=-M/2+1:M/2
        if s^2<=r2
            projection(s+M/2,i)=2*A*B*sqrt(r2-s^2)/r2;
        else
            projection(s+M/2,i)=0;
        end
    end
end

%figure,imshow(projection,[]),title('正弦图');

%直接反投影
bp=zeros(M);
for x=1:M
    for y=1:M
        for i=1:180
            theta=i*pi/180;
            s=(x-M/2)*cos(theta)+(y-M/2)*sin(theta)+M/2;
            si=round(s);
            if si>0&&si<=M
               bp(x,y)=bp(x,y)+projection(si,i);
            end
        end
    end
end
%figure,imshow(bp,[]),title('直接反投影');

subplot(3,3,1),imshow(Img,[]),title('原图')
subplot(3,3,2),imshow(projection,[]),title('原投影正弦图')
subplot(3,3,4),imshow(bp,[]),title('直接反投影')

%滤波反投影
%设置fft的宽度
width=2^nextpow2(size(projection,1));
pro_fft=fft(projection,width);%中心切片定理

%滤波
filter = 2*pi*[0:(width/2-1), width/2:-1:1]'/width;
pro_filter=zeros(width,180);
for i=1:180
    pro_filter(:,i)=pro_fft(:,i).*filter;
end
subplot(3,3,5),imshow(pro_fft,[]),title('傅里叶变换')
subplot(3,3,6),imshow(pro_filter,[]),title('傅里叶变换+滤波')

%ifft求滤波后的投影值
pro_ifft=real(ifft(pro_filter));
subplot(3,3,3),imshow(pro_ifft,[]),title('滤波后的投影正弦图')

%滤波后的反投影
fbp=zeros(M);
for x=1:M
    for y=1:M
        for i=1:180
            theta=i*pi/180;
            s=(x-M/2)*cos(theta)+(y-M/2)*sin(theta)+M/2;
            si=round(s);
%             if si>0&&si<=M
%                 floor=round(s-0.5);
%                 top=round(s+0.5);
%                 if floor>0&&top<=512
%                     a=pro_ifft(floor,i)+(s-floor)*(pro_ifft(top,i)-pro_ifft(floor,i))/(top-floor);
%                     fbp(x,y)=fbp(x,y)+a;
%                 end
%             end
            if si>0&&si<=M
                fbp(x,y)=fbp(x,y)+pro_ifft(si,i);
            end
        end
    end
end
%fbp=fbp/M;
subplot(3,3,7),imshow(fbp,[]),title('滤波反投影图像');


%卷积反投影
%先设计一下|w|在时域中的形式
%采样间隔为1
tau=1;
%采样范围是 -M/2-1:M/2
%p72
g=-(M/2-1):(M/2);

for i=1:M
	if g(i)==0
		hl(i)=1/4*tau^2;
	else if mod(g(i),2)==0
		hl(i)=0;
		else
		hl(i)=-1/(pi^2*tau^2*(g(i)^2));
		end
	end
end

%现在做卷积
k=M/2:3*M/2-1;	%取卷积结果
bp_conv=zeros(M);
for i=1:180
	u=conv(hl,projection(:,i));
	pro_conv=u(k);
	theta=i*pi/180;
	
	for x=1:M
		for y=1:M
			s=(x-M/2)*cos(theta)+(y-M/2)*sin(theta)+M/2;
			s=round(s);
			if s>0&&s<=M
				bp_conv(x,y)=bp_conv(x,y)+pro_conv(s);
			end
		end
	end
end
subplot(3,3,8),imshow(bp_conv,[]),title('卷积反投影')
			
			
			
	


