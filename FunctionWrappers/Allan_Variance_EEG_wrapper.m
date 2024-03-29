function [AllanVar] = Allan_Variance_EEG_wrapper(EEG)

for i=1:EEG.nbchan
    DATA0=EEG.data(i,:)';
    SIZE_DATA0=2^(nextpow2(length(DATA0))-1);
    
    for SCALE=0:17
        S = 2^SCALE;
        
        Z_SIZE = SIZE_DATA0/S;
        
        Z=zeros(Z_SIZE,1);
        
        for id=1:Z_SIZE
            Z(id,1) = mean(DATA0(S*(id-1)+1:S*id));
        end
        
        C = zeros(Z_SIZE-1,1);
        for id=1:Z_SIZE-1
            C(id,1) = (Z(id+1,1) - Z(id,1))^2;
        end
        
        ALLAN_VARIANCE(SCALE+1,i)=mean(C);
    end
    %AllanVar = 
end