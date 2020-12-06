function [dist,ratio,ids,errids, mindists] = dictdist_normc(approxD,D,epsilon)
%DICTDIST Distance between dictionaries.
%  [DIST,RATIO] = DICTDIST(APPROXD,D) computes the distance between the
%  approximate dictionary APPROXD and the true dictionary D, where APPROXD
%  is NxK and D is NxM.
%
%  The distance between the dictionary APPROXD and a single atom A of D is
%  defined as:
%
%      DIST(APPROXD,A) = min  { 1-abs(APPROXD(:,i)' * A) }
%                         i
%
%  The distance between the dictionaries APPROXD and D is defined as:
%
%      DIST(APPROXD,D) = sum { dist(APPROXD, D(:,k)) } / M
%                         k
%
%  Note that 0 <= DIST(APPROXD,D) <= 1, where 0 implies that all atoms in D
%  appear in APPROXD, and 1 implies that the atoms of D are orthogonal to
%  APPROXD.
%
%  The similarity ratio between APPROXD and D is defined as:
%
%      RATIO(APPROXD,D) = #(atoms in D that appear in APPROXD) / M
%
%  where two atoms are considered identical when DIST(A1,A2) < EPSILON with
%  EPSILON=0.01 by default. Note that 0 <= RATIO(APPROXD,D) <= 1, where 0
%  means APPROXD and D have no identical atoms, and 1 means that all atoms
%  of D appear in APPROXD.
%
%  [DIST,RATIO] = DICTDIST(DICT1,DICT2,EPSILON) specifies a different value
%  for EPSILON.
%
%  [DIST,RATIO,IDS] = DICTDIST(...) also returns for each column of D the
%  index of its nearest column in APPROXD, according to the distance
%  measure defined above, DIST(APPROXD,A). Note that IDS is not guaranteed
%  to be a permutaiton, and specifically, may contain repetitions.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  October 2007


if nargin < 3
    epsilon = 0.01;
end
epsilon=max(epsilon,0.01);

[dummy,m] = size(D);

% approxD = normcols(approxD);
% D = normcols(D);

% jing -- replace 'normcols' to 'normc'
approxD = normc(approxD);
D = normc(D);

ids = zeros(1,m);
mindists = zeros(1,m);
for i = 1:m
    atom = D(:,i);
    distances = 1-abs(atom'*approxD);
    [mindists(i),ids(i)] = min(distances);
end
dist = sum(mindists)/m;
cond = mindists<epsilon;
ratio = sum(cond)/m;
errids = find(cond==0);
end