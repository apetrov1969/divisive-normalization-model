function retval_indexes = FindClosestIdx(matrix,target)
% FindClosest -- Find a value closest to 'target' and its index

[dummy,ind] = min(abs(matrix-target));
retval_indexes = ind;

end