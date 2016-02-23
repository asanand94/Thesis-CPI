function A = chungLuP(N,p,q,r)
% Chung-Lu graphs with power law weights
% Here, p=0.1, q=0 ==> Complete graph!
% p (1.0  - 0.9).9
% q (0.0  - 1.0).9
% r (0.01 - 0.5).5

w=(0:N-1);
w=N*p*(1-q*w/N).^(r);
N=size(w,2);
D=w'*w/sum(w);
A = triu(rand(N)<D,1);
A = A+A';