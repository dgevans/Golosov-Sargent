% The following starting values provide a rough solution.
clear all
clc
x = -ones(9, 1);
matlabpool open
tic
parfor i=1:10000
[xOut, fvec, user, ifail] = c05qb(@fcn05qb, x);
end
toc

tic
parfor i=1:10000
[xOut, fvec, ifail] = c05nb(@fcn05nb, x);
end
toc
matlabpool close