function retval_value = ValueFloatIdx(matrix,idx)
% FindClosest -- Find a value closest to 'target' and its index

maxI = size(matrix(:));
maxI = maxI(1,1);

if isnan(idx)
    retval_value = NaN;
elseif idx>maxI
    % retval_value = matrix(maxI);
    retval_value = NaN;
elseif idx<1
    % retval_value = matrix(1);
    retval_value = NaN;
elseif floor(idx)==ceil(idx)
    retval_value = matrix(ceil(idx));
else
    i1=floor(idx);
    i2=ceil(idx);
    w1=i2-idx;
    w2=idx-i1;
    retval_value = w1 * matrix(i1)...
                 + w2 * matrix(i2);
end

end