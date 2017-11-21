function [retval_Max, retval_i1, retval_i2] = Max2D(mtrx)

[list_max,list_i] = max(mtrx);

[retval_Max, retval_i2] = max(list_max);
retval_i1 = list_i(retval_i2);