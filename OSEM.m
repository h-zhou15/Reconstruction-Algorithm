clc,clear all;

%定义系统传输矩阵
CIJ=zeros(128*128,128,60);

%对于每一个角度而言会又一个128*128的矩阵
Nx=129;
Ny=129;

%定义X，Y轴网格的位置
Xplane=[0:1:128]';
Yplane=[0:1:128]';

%定义探测器的初始位置
X10=-36;Y10=0;
X20=164;Y20=0;
D12=200;

%定义初始Alphax,Alphay
Ax=zeros(1,Nx);
Ay=zeros(1,Ny);

for theta =0:6:354
	for j=1:128
		%1，2点随角度的变化关系
		X1=X10*cos(pi*theta /180)-(j-1)*sin(pi*theta/180);
		Y1=X10*sin(pi*theta /180)+(j-1)*cos(pi*theta/180);
		X2=X20*cos(pi*theta /180)-(j-1)*sin(pi*theta/180);
		Y2=X20*sin(pi*theta /180)+(j-1)*cos(pi*theta/180);

		%计算Alphamin和Alphamax
		if X2~=X1
			Ax(1)=(Xplane(1)-X1)/(X2-X1);
			Ax(Nx)=(Xplane(Nx)-X1)/(X2-X1);
		elseif X2==X1
			Ax(:)=0;
		end

		if Y2~=Y1
			Ay(1)=(Yplane(1)-Y1)/(Y2-Y1);
			Ay(Ny)=(Yplane(Ny)-Y1)/(Y2-Y1);
		elseif Y2==Y1
			Ay(:)=0;
		end

		amin=[0,min(Ax(1),Ax(Nx)),min(Ay(1),Ay(Ny))];
		amax=[1,max(Ax(1),Ax(Nx)),max(Ay(1),Ay(Ny))];
		Alphamin=max(amin);
		Alphamax=min(amax);
		
		if Alphamax>Alphamin
			%计算相交线的位置，imin, imax ,jmin,jmax
			if (X2-X1)>=0
				imin=Nx-(Xplane(Nx)-Alphamin*(X2-X1)-X1);
				imax=1+(X1+Alphamax*(X2-X1)-Xplane(1));
			elseif (X2-X1)<=0
				imin=Nx-(Xplane(Nx)-Alphamax*(X2-X1)-X1);
				imax=1+(X1+Alphamax*(X2-X1)-Xplane(1));
			end

			if (Y2-Y1)>=0
				jmin=Ny-(Yplane(Ny)-Alphamin*(Y2-Y1)-Y1);
				jmax=1+(Y1+Alphamax*(Y2-Y1)-Yplane(1));
			elseif (Y2-Y1)<=0
				jmin=Ny-(Yplane(Ny)-Alphamax*(Y2-Y1)-Y1);
				jmax=1+(Y1+Alphamax*(Y2-Y1)-Yplane(1));
			end

			AlphaX=zeros(1,round(imax-imin+0.5));
			%计算Alphax,Alphay
			if (X2-X1)>0
				for ix=1:(imax-imin+1)
					AlphaX(ix)=(Xplane(round(ix+imin-1))-X1)/(X2-X1);
				end
			elseif (X2-X1)<0
				for ix=1:(imax-imin+1)
					AlphaX(ix)=(Xplane(round(imax+1-ix))-X1)/(X2-X1);
				end
			end

			AlphaY=zeros(1,round(jmax-jmin+0.5));
			%计算Alphax,Alphay
			if (Y2-Y1)>0
				for jy=1:(jmax-jmin+1)
					AlphaY(jy)=(Yplane(round(jy+jmin-1))-Y1)/(Y2-Y1);
				end
			elseif (Y2-Y1)<0
				for jy=1:(jmax-jmin+1)
					AlphaY(jy)=(Yplane(round(jmax+1-jy))-Y1)/(Y2-Y1);
				end
			end

			%merge AlphaX,AlphaY 并按照升序合并成一个Alpha数组
			%Alpha里面还应该包含alphamin，和alphamax
			%Alpha=zeros(1,(imax-imin+1)+(jmax-jmin+1)+2);
			Alpha=[Alphamin,AlphaX,AlphaY,Alphamax];
			Alpha=sort(Alpha);	%升序排列
			Alen=size(Alpha,2);

			%现在可以计算相交线长度
			l=zeros(1,Alen);
			for m=2:Alen
				l(m)=D12*(Alpha(m)-Alpha(m-1));
				Alphamid=(Alpha(m)+Alpha(m-1))/2;
				I(m)=1+(X1+Alphamid*(X2-X1)-Xplane(1));
				J(m)=1+(Y1+Alphamid*(Y2-Y1)-Yplane(1));
			end

			%存入传输矩阵
			%第i个像素的位置
			for m=2:Alen
				col=fix(I(m))+1;	%交点的列数
				row=fix(J(m))+1;	%行数
				if row>0&&row<=128&&col>0&&col<=128	
					num=128*(row-1)+col;	%表示第几个像素
					CIJ(num,j,theta/6+1)=l(m);	%第i个像素对第j条投影线的相交线长度
				end
			end
		end
	end
end

%迭代过程

%读入数据
filename='whole_bone.raw';
pid =fopen (filename,'r');
data=fread(pid,128*128*60,'float32');
fclose(pid);
data=reshape(data,[128,128,60]);
f=zeros(128,128,128);

fi=ones(1,128*128);	%给一个最初的猜测值
%分成128张切片
for col=1:128
	for s=1:128
		for i=1:60
			pro(s,i)=data(s,col,i);
		end
	end
	
	prosub=zeros(32,60,4);
	
	%对每一张切片做OSEM	
	%按照投影角度来划分子集，这里分为4个子集，每个4°为一个子集
	
	%subset{1,5,9 ,13,17,21,25,29,33,37,41,45,49,53,57}
	%subset{2,6,10,14,18,22,26,30,34,38,42,46,50,54,58}
	%subset{3,7,11,15,19,23,27,31,35,39,43,47,51,55,59}
	%subset{4,8,12,16,20,24,28,32,36,40,44,48,52,56,60}
	
	
	%做迭代，迭代次数为60次
	
	for iter=1:15
		for sub=1:4
			%构建每个子集的传输矩阵
			%每一条投影线都有一个传输矩阵为128*128（在siddon中是一个128*128个一维向量，存储顺序是
			%逐像素从右到左）
			%然后分为4个子集的话，应该有15*128条投影线，所以每次迭代的传输矩阵应该是（128*128，128*15的二维矩阵
			Cij=zeros(128*128,128*15);
			prosub=zeros(1,128*15);
			%分配Cij的矩阵 和投影
			
			for j=1:15
				for i=1:128
					Cij(:,i+128*(j-1))=CIJ(:,i,sub+4*(j-1));
					prosub(i+128*(j-1))=pro(i,sub+4*(j-1));
				end
			end
			
			%上面已经分好子集了，然后应该按照EMML去做迭代
			p=fi*Cij;
			b=prosub./p;
			bp=Cij*(b');
			c=bp./(sum(Cij,2));
			fi=fi.*(c');
		end
    end
	%迭代结束后，存入矩阵中
	c=reshape(fi,128,128)';
	f(:,:,col)=c;
end

%写入raw文件
result='OSEM-zh.raw';
fid=fopen(result,'w+');
cnt=fwrite(fid,f,'float32');
fclose(fid);