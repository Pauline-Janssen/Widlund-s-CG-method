clear, clc, clear all;
N = 50; %size of matrix
alpha = 0.9; 
tol = 0.01; % tolerance for singular values
x = ones(N,1); % given solution of linear system for numerical experiments
limit = 1000; % maximum of iterations
reorth = 0; % choose 1 for double reorthogonalisation, 0 without

L = diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1) + alpha*(diag(1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1));
%L = diag(3*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1) + alpha*(diag(1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1));

% different functions for the numerical experiments
CGW_delay_of_conv(L, x, alpha, N, limit, reorth);
CGW_loss_of_orth_final(L, x, N, limit, reorth)
CGW_rank_defi_final(L, x, alpha, N, limit, reorth, tol)


clear;
N = 50;
alpha = 0.9;
tol = 0.01;
x = ones(N,1); 
limit = 1000;
reorth = 1;

L = diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1) + alpha*(diag(1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1));
%L = diag(3*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1) + alpha*(diag(1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1));

% different functions for the numerical experiments
CGW_delay_of_conv(L, x, alpha, N, limit, reorth);
CGW_loss_of_orth_final(L, x, N, limit, reorth)
CGW_rank_defi_final(L, x, alpha, N, limit, reorth, tol)

%----------------------------------------------------------
%for loop for certain alpha, Chapter 4.1 Influence of the skew-symmetric part
% alphas = [0.01, 0.5, 0.9, 3, 50] ;
% for alpha = 1:length(alphas)  
      % L = diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1) + alphas(alpha)*(diag(1*ones(1,N-1),1) + diag(-1*ones(1,N-1),-1));
      % line = alpha;
      % CGW_scatter(L,line,alphas(alpha),N)
      % CGW_alphas(L, x, alpha, N, tol, limit); 
% end



%----------------------------------------------------------
%  parameter for the experiments in Chapter 4.4 Sharpness of error bounds
%  gamma = 1;
%  skew = build_skew_matrix(N, 100, gamma);
%  L = eye(N)-skew;
%  line = 1,
%  CGW_error_bound(L, x, N, limit, gamma, reorth);
%  CGW_scatter(L, line, gamma, N)

%  gamma = 0.5;
%  skew = build_skew_matrix(N, 100, gamma);
%  L = eye(N)-skew;
%  line = 2;
%  CGW_error_bound(L, x, N, limit, gamma, reorth);
%  CGW_scatter(L, line, gamma, N)


% for loop for matrices with gamma in Chapter 4.4 Sharpness of error bounds
% gammas = [0.01, 0.5, 0.8, 0.9, 1] ;
% for gamma = 1:length(gammas)
%     line = gamma;
%     skew = build_skew_matrix(N,100,gammas(gamma));
%     L = eye(N)-skew; % 10*skew
%     CGW_error_bound(L, x, N, limit, gammas(gamma), reorth);
%     CGW_scatter(L,line,gammas(gamma),N)
% end
