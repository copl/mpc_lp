
m = 10;
n = 100;
A = sprandn(m,n,0.8);
b = A*rand(n,1);
c = A'*randn(m,1) + rand(n,1);

opts = struct; %Empty opts
opts.verbose = true;
opts.max_iter        = 100 ;
opts.ini_mehrotra    = true;
opts.secord          = true;
opts.verbose         = true; 

[x,y,s,info] = mehrotra_lp_solver(A,b,c,opts);

