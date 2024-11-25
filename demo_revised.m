%miu
%theta
clear;close all;clc
I = [11 15 17 500]; % tensor dimension
R = [3,5,2,4];      % TR rank 
R = [R R(1)];
N = length(I);
factor_noise = 0e-3;
Z_true = cell(1,3);
for k = 1 : N
    Z_true{k} = randn(R(k),I(k),R(k+1));
end
% Generate Tensor T From TR-Cores
Y = coreten2tr(Z_true) + factor_noise * randn(I);
missing_ratio = 0.0;
P_Omega = rand(size(Y));
P_Omega = 1.*(P_Omega <= 1-missing_ratio);

%% Algorithms
miu   = 0.1;
theta = 1e-3;
R = R(2:end);
[G, error] = APTR_SGD(Y, P_Omega,miu, theta,R,Z_true);
figure; plot(error,'b','LineWidth',2); 
ylabel('ERROR'); xlabel('Iteration (aka time index)');
title('Convergence rate')
set(gca, 'YScale', 'log');
 
