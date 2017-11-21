function retval_index = FindFloatIdx(vector,target)
% FindClosest -- Find a value closest to 'target' and its index

vector = vector(:);
maxI = size(vector);
maxI = maxI(1,1); % Vector to Scalar

if target >= max(vector(:))
    retval_index = NaN;
elseif target <= min(vector(:))
    retval_index = NaN;
else
    for i=1:maxI-1
        if ( (vector(i)-target)*(vector(i+1)-target) <= 0 )
            break;
        end
    end

    if vector(i) == vector(i+1)
        retval_index = i+0.5;
    else
        retval_index = (i)  *( vector(i+1)-target)/(vector(i+1)-vector(i))...
                      +(i+1)*(-vector(i)  +target)/(vector(i+1)-vector(i));
    end
end

end