function [v]=global_vel(a,b,e,f,g,v,gamma,n,alpha,U)
for j=1:n-1
        [P]=influence_matrix_vortex(a(j),a(j+1),b(j),b(j+1),e(j),f,g);
        v=v+P*[gamma(j),gamma(j+1)]';
end  
v=[U*cos(alpha),U*sin(alpha)]'+v;
