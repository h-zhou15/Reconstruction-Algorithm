clc,clear all;

%读入数据
filename='whole_bone.raw';
pid =fopen (filename,'r');
data=fread(pid,128*128*60,'float32');
fclose(pid);
data=reshape(data,[128,128,60]);

%把投影数据沿着纸面向里（即就是沿着探测器的列方向做fbp）的方向分成128张切片

M=128;

fbp=zeros(128,128,128);
pro=zeros(128,60);
pro_ifft=zeros(128,60,128);

%分成128张切片
for col=1:128
	for s=1:128
		for i=1:60
			pro(s,i)=data(s,col,i);%每一张切片的投影信息
		end
	end
	
	%下面是对一张切片做fbp
	%基本思路：fft->滤波->ifft->反投影
	
	width=size(pro,1);	%设置fft的宽度
	for i=1:60
		pro_fft(:,i)=fft(pro(:,i),width);
	end
	
	%滤波,选用斜坡滤波
	filter=2*[0:width/2-1,width/2:-1:1]'/width;	%|w|
	for i=1:60
		pro_filter(:,i)=pro_fft(:,i).*filter;
	end
	
	%逆傅里叶变换得到滤波后的投影值
	pro_ifft(:,:,col)=real(ifft(pro_filter));
	
	%反投影
	for x=1:128
		for y=1:128
			for i=1:60
				theta=6*(i-1)*pi/180;
				s=M/2+(x-M/2)*cos(theta)+(y-M/2)*sin(theta);
				%考虑到反投影后的s可能落在两个探测器中间
				%则采用线性插值的办法计算得到每次累加的投影值
				flo=round(s-0.5);
				top=round(s+0.5);
				
				if top<=128&&flo>0&&s~=fix(s)
					p=pro_ifft(flo,i,col)+(s-flo)*(pro_ifft(top,i,col)-pro_ifft(flo,i,col));%线性插值
                elseif s==fix(s)&&s>0&&s<=128
					p=pro_ifft(s,i,col);
				end
				
				%反投影
				if s>0&&s<=128
					fbp(x,y,128-col+1)=fbp(x,y,128-col+1)+p;    	
				end
			end
		end
	end
end


%写入raw文件
result='FBP-zh.raw';
fid=fopen(result,'w+');
cnt=fwrite(fid,fbp,'float32');
fclose(fid);

%图像评价
%data为原始投影值
%pro_ifft为滤波后的投影值
RDev=zeros(128);
%找到每一张切片的最大值

for col=1:128
	for x=1:128
		for i=1:60
			fmax=max(max(data(:,col,:)));
			RDev(col)=RDev(col)+((data(x,col,i)-pro_ifft(x,i,128-col+1))/fmax)^2;
		end
	end
end
fprintf('mean of RDev : %f\n',mean(RDev(:,1)));
fprintf('the error for single pixel: %f \n',mean(RDev(:,1))/(128*128));

