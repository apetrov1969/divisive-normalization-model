function [retval_values,retval_indexes] = FindClosest(matrix,target)
% FindClosest -- Find a value closest to 'target' and its index

[dummy,ind] = min(abs(matrix-target));

retval_values  = matrix(ind);
retval_indexes = ind;

end