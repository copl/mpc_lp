clear all
close all

%This is a Mehrotra predictor-corrector implementation for LP

n = 100; %Size of the cone
m = 50; %Number of costraints
nu = n; %Complexity of the barrier

%-------
%Generate a feasible primal point
x = rand(n,1);

%Generate a feasible dual point
s = rand(n,1);

%Generate a feasible and bounded primal dual problem
A = sprandn(m,n,0.5);
b = A*x;
s = rand(n,1);
c = A'*ones(m,1) + s;
fprintf('Generated a random feasible LP with %i constraints and %i variables\n',m,n);
clear x s

opts = struct; %Empty opts
opts.verbose = true;
opts.max_iter        = 100 ;
opts.ini_mehrotra    = true;
opts.secord          = true;
opts.verbose         = true;
[x,y,s,info] = mehrotra_lp_solver(A,b,c,opts);
