function [X,err] = AnalysisL1Denoising(Y,Omega,lambda)
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
X = Y;
Z = Omega*X;
gamma = lambda; % ALMM Lagrange multiplier
i = 1;
Theta = zeros(N,L); % ALMM Lagrange multipliers
LM = inv(lambda*eye(M)+ gamma*(Omega')*Omega);
while ((i <= 5) || (norm(Z-Omega*X,'fro') >= .001)) && (i <= 1000)
    
%     i
    
    X = LM*(lambda*Y + gamma*Omega'*(Z-Theta));
    Z = softthreshold(Omega*X+Theta,1/gamma);
    err(i) =  sum(sum(abs(Z)))+lambda/2 *(norm(X-Y,'fro')^2);
    Theta = Theta - (Z-Omega*X); 
    i = i+1;
end