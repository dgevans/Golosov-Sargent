% The following starting values provide a rough solution.
clear all
clc
x = -ones(9, 1);
tic
system('module load nag/mbl6a22dml')
for i=1:1000
[xOut, fvec, ifail] = c05nb(@fcn05nb, x);
end
toc

system('module unload nag')
system('module load nag/mbl6a23dml')
tic
for i=1:1000
[xOut, fvec, user, ifail] = c05qb(@fcn05qb, x);
end
toc

for i=1:1000
[xOut, fvec, ifail] = c05nb(@fcn05nb, x);
end
toc
