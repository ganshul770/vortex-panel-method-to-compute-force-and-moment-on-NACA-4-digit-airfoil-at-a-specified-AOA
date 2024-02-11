clear all;
close all;
clc;

%%
%naca 2412 nomenclature
u=input('first digit of naca = ');
v=input('second digit of naca = ');
w=input('last two digit of naca = ');
ymc_c=u*0.01;
xmc_c=v*0.1;
tm=w*0.01;
c=1;
z=input("alpha value = ");
alpha=z*pi/180;
U=10;

%%
%cosine clustring
%choose even number of points on circle 
n=160;
theta=2*pi/(n-1);
i=1:n/2;
x_c=0.5*(1-cos((i-0.5)*theta));
[xl,yl,xu,yu]=surface_coordinate(tm,xmc_c,ymc_c,c,x_c);
x=[fliplr(xl),xu];
y=[fliplr(yl),yu];
for i=1:n-1
    dx(i)=x(i+1)-x(i);
    dy(i)=y(i+1)-y(i);
    x_cp(i)=(x(i+1)+x(i))/2;
    y_cp(i)=(y(i+1)+y(i))/2;
    l(i)=sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end

%%
A=zeros(n,n);
T=zeros(n,n);
B=zeros(n,1);
v_in_t=zeros(n,1);
for i=1:n-1
    for j=1:n-1
        [P]=influence_matrix_vortex(x(j),x(j+1),y(j),y(j+1),l(j),x_cp(i),y_cp(i));
        A(i,j)=A(i,j)+(-P(1,1)*(y(i+1)-y(i))/l(i)+P(2,1)*(x(i+1)-x(i))/l(i));
        A(i,j+1)=A(i,j+1)+(-P(1,2)*(y(i+1)-y(i))/l(i)+P(2,2)*(x(i+1)-x(i))/l(i));
    end  
end
for i=1:n-1
    for j=1:n-1
        [P]=influence_matrix_vortex(x(j),x(j+1),y(j),y(j+1),l(j),x_cp(i)+0.001*(-(y(i+1)-y(i))/l(i)),y_cp(i)+0.001*((x(i+1)-x(i))/l(i)));
        T(i,j)=T(i,j)+(P(1,1)*(x(i+1)-x(i))/l(i)+P(2,1)*(y(i+1)-y(i))/l(i));
        T(i,j+1)=T(i,j+1)+(P(1,2)*(x(i+1)-x(i))/l(i)+P(2,2)*(y(i+1)-y(i))/l(i));
    end  
end

A(n,[1,end])=1;
T(n,[1,end])=1;
for i=1:n-1
    B(i,1)=U*((y(i+1)-y(i))*cos(alpha)-(x(i+1)-x(i))*sin(alpha))/l(i);
    v_in_t(i,1)=U*((x(i+1)-x(i))*cos(alpha)+(y(i+1)-y(i))*sin(alpha))/l(i);
end
gamma=A\B;
T_vel=T*gamma+v_in_t;

%%
%cp variation
cp=1-(T_vel.^2./U^2);
figure(1)
plot(x_cp,cp(1:n-1));
xlabel('---- x/c ---->');
ylabel('---- cp ---->');
title("cp on upper and lower surface");
set(gca, 'YDir','reverse');

%% 
% naca profile
figure(2)
plot(x,y);
xlabel('---- x/c ---->');
ylabel('---- y/c ---->');
title('naca 2412');

%%
cl=0;
cm_le=0;
for j=1:n-1
    cl=cl+l(j)*(gamma(j)+gamma(j+1))/(U*c);
    cm_le=cm_le-l(j)*((2*x(j)*gamma(j)+x(j)*gamma(j+1)+x(j+1)*gamma(j)+2*x(j+1)*gamma(j+1))*cos(alpha)+(2*y(j)*gamma(j)+y(j)*gamma(j+1)+y(j+1)*gamma(j)+2*y(j+1)*gamma(j+1))*sin(alpha))/(3*U*c^2);
end
cm_qc=cm_le+cl*(0.25);


%%
% streamline
%from trailing edge
% t_T(1)=1;
% s_T(1)=0;
% k=[1;0];
% for i=2:10
%     v=zeros(2,1);
%     [v]=global_vel(x,y,l,t_T(i-1),s_T(i-1),v,gamma,n,alpha,U);
%     %euler method
%     k=k+0.01*v;
%     t_T(i)=k(1);
%     s_T(i)=k(2);
% end
    [t_T,s_T]=sl_by_point(x,y,l,gamma,n,alpha,U,1,0,0.01,30);


%%
a=[x(0.5*n-1),0,x(0.5*n+1)];
b=[y(0.5*n-1),0,y(0.5*n+1)];
vt=zeros(3,1);
for i=1:3
    h=zeros(2,1);
    [h]=global_vel(x,y,l,a(i)+10^(-4),b(i)+10^(-4),h,gamma,n,alpha,U);

    vt(i)=[(x(0.5*n-1+i)-a(i))/sqrt((x(0.5*n-1+i)-a(i))^2+(y(0.5*n-1+i)-b(i))^2),(y(0.5*n-1+i)-b(i))/sqrt((x(0.5*n-1+i)-a(i))^2+(y(0.5*n-1+i)-b(i))^2)]*h;
end
d_u_vt=-(vt(3)-vt(2))/(a(3)-a(2));
d_l_vt=(vt(2)-vt(1))/(a(2)-a(1));


%%
vt_le=vt(2);
N_R=0;
u_vt_n=0;
l_vt_n=0;
while (u_vt_n<10^(-5)*U || l_vt_n<10^(-5)*U)
    if (vt_le<0)%search upper
        N_R=N_R-vt_le/d_u_vt;
    else%search lower
        N_R=N_R-vt_le/d_l_vt;
    end
    xl_g=0;yl_g=0;xu_g=0;yu_g=0;
    [xl_n,yl_n,xu_n,yu_n]=surface_coordinate(tm,xmc_c,ymc_c,c,N_R);
    h_u=zeros(2,1);
    h_l=zeros(2,1);
    [h_u]=global_vel(x,y,l,xu_n,yu_n,h_u,gamma,n,alpha,U);
    [h_l]=global_vel(x,y,l,xl_n,yl_n,h_l,gamma,n,alpha,U);
    xl_n=[xl_g,xl_n];yl_n=[yl_g,yl_n];xu_n=[xu_g,xu_n];yu_n=[yu_g,yu_n];
    u_vt_n=[(xu_n(2)-xu_n(1))/sqrt((xu_n(2)-xu_n(1))^2+(yu_n(2)-yu_n(1))^2),(yu_n(2)-yu_n(1))/sqrt((xu_n(2)-xu_n(1))^2+(yu_n(2)-yu_n(1))^2)]*h_u;
    l_vt_n=[(xl_n(2)-xl_n(1))/sqrt((xl_n(2)-xl_n(1))^2+(yl_n(2)-yl_n(1))^2),(yl_n(2)-yl_n(1))/sqrt((xl_n(2)-xl_n(1))^2+(yl_n(2)-yl_n(1))^2)]*h_l;
    xl_g=xl_n(2);yl_g=yl_n(2);xu_g=xu_n(2);yu_g=yu_n(2);
end

if u_vt_n<l_vt_n
    xs=xu_g;
    ys=yu_g;
    v_s=h_u;
else
    xs=xl_g;
    ys=yl_g;
    v_s=h_l;
end

%%
%from leading edge
[t_L,s_L]=sl_by_point(x,y,l,gamma,n,alpha,U,xs,ys,-0.001,500);

figure(3)
plot(x,y);
xlabel('---- x/c ---->');
ylabel('---- y/c ---->');
title("streamline");
hold on;
plot(t_T,s_T);
hold on;
plot(t_L,s_L);
hold on;
%other streamlines
for i=-29:2:29
[m,o]=sl_by_point(x,y,l,gamma,n,alpha,U,t_L(end),s_L(end)+0.03*i,0.01,100);
plot(m,o);
hold on;
end
