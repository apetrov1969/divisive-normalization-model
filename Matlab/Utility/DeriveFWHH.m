function [retval_fullH,retval_fullH_idx,...
          retval_halfH,retval_halfH_idx] = DeriveFWHH(curve)

retval_fullH = NaN;
retval_halfH = NaN;
retval_fullH_idx = NaN;
retval_halfH_idx = NaN;

curve = curve(:);

[fullH,fullH_idx] = max(curve(:));

if(size(fullH_idx(:),1)~=1)
    return;
end

if(fullH_idx==1 || fullH_idx==size(curve,1))
    return;
end

curveL = curve(1:fullH_idx-1);
curveR = curve(fullH_idx+1:size(curve,1));
halfH = fullH/2.0;

%[halfHL,halfHL_idx] = FindClosest(curveL,halfH);
%[halfHR,halfHR_idx] = FindClosest(curveR,halfH);

halfHL_idx = FindFloatIdx(curveL,halfH);
halfHR_idx = FindFloatIdx(curveR,halfH) + fullH_idx;

retval_fullH     = fullH;
retval_halfH     = halfH;
retval_fullH_idx = fullH_idx;
retval_halfH_idx = [halfHL_idx,halfHR_idx];

end %%% of file