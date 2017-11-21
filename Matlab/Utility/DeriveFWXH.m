function [retval_fullWhalfH,...
          retval_fullH,retval_fullH_idx,...
          retval_xxxxH,retval_xxxxH_idx] = DeriveFWXH(curve,height)

retval_fullWhalfH = NaN;
retval_fullH = NaN;
retval_xxxxH = NaN;
retval_fullH_idx = NaN;
retval_xxxxH_idx = NaN;

[fullH,fullH_idx] = max(curve(:));

if(size(fullH_idx(:),1)~=1)
    return;
end

if(fullH_idx==1 || fullH_idx==size(curve,2))
    return;
end

curveL = curve(1:fullH_idx-1);
curveR = curve(fullH_idx+1:size(curve,2));
xxxxH = fullH*height;
[xxxxHL,xxxxHL_idx] = FindClosest(curveL,xxxxH);
[xxxxHR,xxxxHR_idx] = FindClosest(curveR,xxxxH);
    
if(xxxxHL_idx==1 && xxxxHL_idx==size(curveL,2))
    return;
end

if(xxxxHR_idx==1 || xxxxHR_idx==size(curveR,2))
    return;
end

xxxxHR_idx = xxxxHR_idx+fullH_idx;

retval_fullH     = fullH;
retval_fullH_idx = fullH_idx;
retval_xxxxH     = [xxxxHL,xxxxHR];
retval_xxxxH_idx = [xxxxHL_idx,xxxxHR_idx];
retval_fullWhalfH = xxxxHR_idx-xxxxHL_idx;

end %%% of file