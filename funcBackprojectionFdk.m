function [rec] = funcBackprojectionFdk(Q,pixelsize,SOD,beta,beta_num,N,detector_size,SDD)
    rec=zeros(N,N,N);
    t=SDD/SOD;
    N_d=size(Q,1);
    for k1=1:N
        x=pixelsize*(k1-(N+1)/2);
        for k2=1:N
            y=pixelsize*((N+1)/2-k2);
            U=(SOD+y*sin(beta)+x*cos(beta))/SOD;
            a=t*(y*cos(beta)-x*sin(beta))/U;
         %   X=round(a/detector_size+N_d/2+1);
            xx=floor(a/detector_size+N_d/2+1); %整数部分
            u1=a/detector_size+N_d/2+1-xx; %小数部分
            
            for k3=1:N
                z=pixelsize*(k3-(N+1)/2);
                b=t*z/U;
      %          Y=round(b/detector_size+N_d/2+1);
                yy=floor(b/detector_size+N_d/2+1);
                u2=b/detector_size+N_d/2+1-yy;
                
                %双线性插值
                if (xx>=1)&&(xx<N)&&(yy>=1)&&(yy<N)
 %              if (X>=1&&X<=N)&&(Y>=1&&Y<=N)
                    temp=(1-u1)*(1-u2)*Q(xx,yy)+(1-u1)*u2*Q(xx,yy+1)+u1*(1-u2)*Q(xx+1,yy)+u1*u2*Q(xx+1,yy+1);
 %                   temp=Q(X,Y);
                    rec(k2,N+1-k3,k1)=rec(k2,N+1-k3,k1)+temp/(U.^2)*(2*pi/beta_num);
                else
                    rec(k2,N+1-k3,k1)=rec(k2,N+1-k3,k1);
                end
            end
        end
    end
    
end

