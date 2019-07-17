%step 3
%2D diffusion equationbs
% Note: (mu*dt)/dx^2 has to be <1


clear all;
close all;
clc;

%constant
mu=0.3;

%Mesh
xl=2;yl=2;
x1=0.5;y1=0.5;
x2=1;y2=1;
dx=0.025;dy=0.025;
Nx=round(xl/dx);
Ny=round(yl/dy);
x=0:dx:xl;
y=0:dy:yl;

%Initial condition
dt=0.0005;
tf=0.1;
nt=round(tf/dt);
size_x=size(x);size_y=size(y);
u=ones(size_x(2),size_y(2));
nx1=round(x1/dx);ny1=round(y1/dy);
nx2=round(x2/dx);ny2=round(y2/dy);
u(nx1:nx2,ny1:ny2)=2;

%loop

for t=1:nt
    for i=2:Nx-1
        for j=2:Ny-1
        u(i,j)=u(i,j)+((mu*dt)/dx.^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))+((mu*dt)/dy.^2)*(u(i,j+1)-2*u(i,j)+u(i,j-1));
        end
    end
      figure(1);
      shading interp;
      axis([0 2 0 2 1 2]);
      pause(0.1);
%       hold on;
      surf(x,y,u);
end
 
