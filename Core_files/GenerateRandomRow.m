function [w] = GenerateRandomRow(x)
%GenerateRandomRow generates one analysis atom from a given data set of 
%signals, by drawing at random d-1 signals and computing the vector that 
%spans their nullspace.
%  [w] = GenerateRandomRow(x)
%  ========================================================================
%  Input:
%  x - matrix containing a set of signals (in its columns).
%  ========================================================================
%  Output:
%  w - row corresponding to one analysis atom (normalized).
%  ========================================================================
%  Tomer Peleg
%  Department of Electrical Engineering
%  Technion, Haifa 32000 Israel
%  tomerfa@tx.technion.ac.il
%
%  October 2012
%  ========================================================================
[d,N]=size(x);
ids1 = randperm(N);
ids1 = ids1(1:d-1);
SubS = x(:,ids1);
w=randn(1,d);
w=w-w*SubS*pinv(SubS);
w=w/sqrt(w*w');