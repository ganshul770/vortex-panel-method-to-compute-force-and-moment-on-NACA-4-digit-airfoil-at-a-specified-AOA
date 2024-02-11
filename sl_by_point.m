function [t_T,s_T]=sl_by_point(x,y,l,gamma,n,alpha,U,a,b,dx,N)
t_T(1)=a;
s_T(1)=b;
% v=zeros(2,1);
k=[a;b];
for i=2:N
    v=zeros(2,1);
    [v]=global_vel(x,y,l,t_T(i-1),s_T(i-1),v,gamma,n,alpha,U);
    %euler method
    k=k+dx*v;
    t_T(i)=k(1);
    s_T(i)=k(2);
end