clc,clear all;

%����ϵͳ�������
CIJ=zeros(128*128,128,60);

%����ÿһ���Ƕȶ��Ի���һ��128*128�ľ���
Nx=129;
Ny=129;

%����X��Y�������λ��
Xplane=[0:1:128]';
Yplane=[0:1:128]';

%����̽�����ĳ�ʼλ��
X10=-36;Y10=0;
X20=164;Y20=0;
D12=200;

%�����ʼAlphax,Alphay
Ax=zeros(1,Nx);
Ay=zeros(1,Ny);

for theta =0:6:354
	for j=1:128
		%1��2����Ƕȵı仯��ϵ
		X1=X10*cos(pi*theta /180)-(j-1)*sin(pi*theta/180);
		Y1=X10*sin(pi*theta /180)+(j-1)*cos(pi*theta/180);
		X2=X20*cos(pi*theta /180)-(j-1)*sin(pi*theta/180);
		Y2=X20*sin(pi*theta /180)+(j-1)*cos(pi*theta/180);

		%����Alphamin��Alphamax
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
			%�����ཻ�ߵ�λ�ã�imin, imax ,jmin,jmax
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
			%����Alphax,Alphay
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
			%����Alphax,Alphay
			if (Y2-Y1)>0
				for jy=1:(jmax-jmin+1)
					AlphaY(jy)=(Yplane(round(jy+jmin-1))-Y1)/(Y2-Y1);
				end
			elseif (Y2-Y1)<0
				for jy=1:(jmax-jmin+1)
					AlphaY(jy)=(Yplane(round(jmax+1-jy))-Y1)/(Y2-Y1);
				end
			end

			%merge AlphaX,AlphaY ����������ϲ���һ��Alpha����
			%Alpha���滹Ӧ�ð���alphamin����alphamax
			%Alpha=zeros(1,(imax-imin+1)+(jmax-jmin+1)+2);
			Alpha=[Alphamin,AlphaX,AlphaY,Alphamax];
			Alpha=sort(Alpha);	%��������
			Alen=size(Alpha,2);

			%���ڿ��Լ����ཻ�߳���
			l=zeros(1,Alen);
			for m=2:Alen
				l(m)=D12*(Alpha(m)-Alpha(m-1));
				Alphamid=(Alpha(m)+Alpha(m-1))/2;
				I(m)=1+(X1+Alphamid*(X2-X1)-Xplane(1));
				J(m)=1+(Y1+Alphamid*(Y2-Y1)-Yplane(1));
			end

			%���봫�����
			%��i�����ص�λ��
			for m=2:Alen
				col=fix(I(m))+1;	%���������
				row=fix(J(m))+1;	%����
				if row>0&&row<=128&&col>0&&col<=128	
					num=128*(row-1)+col;	%��ʾ�ڼ�������
					CIJ(num,j,theta/6+1)=l(m);	%��i�����ضԵ�j��ͶӰ�ߵ��ཻ�߳���
				end
			end
		end
	end
end

%��������

%��������
filename='whole_bone.raw';
pid =fopen (filename,'r');
data=fread(pid,128*128*60,'float32');
fclose(pid);
data=reshape(data,[128,128,60]);
f=zeros(128,128,128);

fi=ones(1,128*128);	%��һ������Ĳ²�ֵ
%�ֳ�128����Ƭ
for col=1:128
	for s=1:128
		for i=1:60
			pro(s,i)=data(s,col,i);
		end
	end
	
	prosub=zeros(32,60,4);
	
	%��ÿһ����Ƭ��OSEM	
	%����ͶӰ�Ƕ��������Ӽ��������Ϊ4���Ӽ���ÿ��4��Ϊһ���Ӽ�
	
	%subset{1,5,9 ,13,17,21,25,29,33,37,41,45,49,53,57}
	%subset{2,6,10,14,18,22,26,30,34,38,42,46,50,54,58}
	%subset{3,7,11,15,19,23,27,31,35,39,43,47,51,55,59}
	%subset{4,8,12,16,20,24,28,32,36,40,44,48,52,56,60}
	
	
	%����������������Ϊ60��
	
	for iter=1:15
		for sub=1:4
			%����ÿ���Ӽ��Ĵ������
			%ÿһ��ͶӰ�߶���һ���������Ϊ128*128����siddon����һ��128*128��һά�������洢˳����
			%�����ش��ҵ���
			%Ȼ���Ϊ4���Ӽ��Ļ���Ӧ����15*128��ͶӰ�ߣ�����ÿ�ε����Ĵ������Ӧ���ǣ�128*128��128*15�Ķ�ά����
			Cij=zeros(128*128,128*15);
			prosub=zeros(1,128*15);
			%����Cij�ľ��� ��ͶӰ
			
			for j=1:15
				for i=1:128
					Cij(:,i+128*(j-1))=CIJ(:,i,sub+4*(j-1));
					prosub(i+128*(j-1))=pro(i,sub+4*(j-1));
				end
			end
			
			%�����Ѿ��ֺ��Ӽ��ˣ�Ȼ��Ӧ�ð���EMMLȥ������
			p=fi*Cij;
			b=prosub./p;
			bp=Cij*(b');
			c=bp./(sum(Cij,2));
			fi=fi.*(c');
		end
    end
	%���������󣬴��������
	c=reshape(fi,128,128)';
	f(:,:,col)=c;
end

%д��raw�ļ�
result='OSEM-zh.raw';
fid=fopen(result,'w+');
cnt=fwrite(fid,f,'float32');
fclose(fid);