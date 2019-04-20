function [Q] = funcFilter(RF1,fh_RL,N_d)
%沿着a方向做卷积滤波
    for i=1:N_d
        RF2=RF1(:,i);
        k=N_d:2*N_d-1;
        q=conv(RF2,fh_RL);
        Q(:,i)=q(k);
    end
end

