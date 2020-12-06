function [X,cost] = ADxh108_reliable(Y,Omega,reliable_samples_mat,X_sparse,ClippingLevel,clipped_pos_mat,clipped_neg_mat,ZZ)
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

% Z = Omega*Y;
Z=X;
lambda=1;
gamma = lambda; % ALMM Lagrange multiplier
i = 1;
d=size(Y,1);
LM = inv(lambda*eye(M)+ gamma*(Omega')*Omega);%20200913
clipped_pos_mat=clipped_pos_mat;
clipped_neg_mat=clipped_neg_mat;
% Theta = zeros(N,L); 
while(i <= 1) 
% while ((i <= 1) || (norm(Z-X_sparse+Omega*X,'fro') >= .1)) && (i <= 1000)
tt=1;

  X = LM*(lambda*Y+Omega'*X_sparse);
%   X=X-0.25*(-2*Omega'*X_sparse+2*Omega'*Omega*X);


     ResidualMat = ZZ-X ;
     ResidualMat(clipped_pos_mat) = max(ResidualMat(clipped_pos_mat),0);
     ResidualMat(clipped_neg_mat) = min(ResidualMat(clipped_neg_mat),0);

     X(clipped_pos_mat) = X(clipped_pos_mat)+ResidualMat(clipped_pos_mat);
     X(clipped_neg_mat) = X(clipped_neg_mat)+ResidualMat(clipped_neg_mat);     
     

   X(clipped_pos_mat) = min(X(clipped_pos_mat),1) ;
   X(clipped_neg_mat) = max(X(clipped_neg_mat),-1) ;  

    X(reliable_samples_mat)=ZZ(reliable_samples_mat);
  ResidualMat_temp=Z-X;
 
i=i+1;
%   norm(Z-X_sparse+Omega*X,'fro') 
end
%     X(clipped_pos_mat) = min(X(clipped_pos_mat),1) ;
%    X(clipped_neg_mat) = max(X(clipped_neg_mat),-1) ;  
%    Y(clipped_pos_mat) = max(X(clipped_pos_mat),ClippingLevel);
%    Y(clipped_neg_mat) = min(X(clipped_neg_mat),-ClippingLevel);
%    Y(clipped_pos_mat) = X(clipped_pos_mat);
%    Y(clipped_neg_mat) = X(clipped_neg_mat);

%    X=Y;
% cost = sum(sum(ResidualMat_temp.^2)); 
cost = sum(sum(ResidualMat.^2)); 

 

 
 