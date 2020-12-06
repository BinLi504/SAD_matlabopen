function [X,cost] = SAD_gradient(Y,Omega,reliable_samples_mat,X_sparse,ClippingLevel,clipped_pos_mat,clipped_neg_mat,ZZ)
% Denoising using L1 analysis model, given the operator
% Y = = Noisy Cosparse Signals
% Omega = Analysis Operator
% lambda = Lagrange Multiplier of the Objective
% X = Denoised Signals
% err = The Cost Function Value in Different Iterations

% Constrained Analysis Operator Learning for Cosparse Signal Modelling
% Written by Mehrdad Yaghoobi, version 1.0                                    
%     
% Copyright 2013 Mehrdad Yaghoobi, Sangnam Nam, Remi Gribonval and Mike E. Davies
% 
% For all details please refer to README.m
%
% This software is a free software distributed under the terms of the GNU 
% Public License version 3 (http://www.gnu.org/licenses/gpl.txt). You can 
% redistribute it and/or modify it under the terms of this licence, for 
% personal and non-commercial use and research purpose. 

[N,M] = size(Omega);
[M,L] = size(Y);
X=Y;
Q=Y;

Z = Omega*Y;

lambda=1;
gamma = lambda; % ALMM Lagrange multiplier
i = 1;
d=size(Y,1);

clipped_pos_mat=clipped_pos_mat;
clipped_neg_mat=clipped_neg_mat;
Theta = zeros(N,L); 
while(i <= 10) 

%   X = LM*(lambda*Y+Omega'*X_sparse);
X=X-0.25*(-2*Omega'*X_sparse+2*Omega'*Omega*X);
   ResidualMat = ZZ-X ; 
   ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
   ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);

   X(clipped_pos_mat) = X(clipped_pos_mat)+ResidualMat(clipped_pos_mat);
   X(clipped_neg_mat) = X(clipped_neg_mat)+ResidualMat(clipped_neg_mat);     
     

   X(clipped_pos_mat) = min(X(clipped_pos_mat),1) ;
   X(clipped_neg_mat) = max(X(clipped_neg_mat),-1) ;  
   X(reliable_samples_mat)=ZZ(reliable_samples_mat);

i=i+1;

end

cost = sum(sum(ResidualMat.^2)); 

 

 
 