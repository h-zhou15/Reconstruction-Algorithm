function [RF1] = funcWeightProjectData(RF,N_d,SOD,dd)
    for i=1:N_d
        for j=1:N_d
            RF1(i,j)=RF(i,j)*SOD/sqrt((SOD)^2+((i-(N_d)/2)*dd).^2+(dd*(j-(N_d)/2)).^2);
        end
    end 
end

