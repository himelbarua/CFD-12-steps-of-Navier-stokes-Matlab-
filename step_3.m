%step 3
%1D diffusion equationbs
% Note: (mu*dt)/dx^2 has to be <1


clear all;
close all;
clc;

%constant
mu=0.3;

%Mesh
xl=2;
x1=0.5;
x2=1;
dx=0.1;
N=round(xl/dx);
x=0:dx:xl;

%Initial condition
dt=0.01;
tf=1;
nt=round(tf/dt);
size_x=size(x);
u=ones(size_x(2),1);
n1=round(x1/dx);
n2=round(x2/dx);
u(n1:n2)=2;

%loop

for t=1:nt
    for i=2:N-1
        u(i)=u(i)+((mu*dt)/dx.^2)*(u(i+1)-2*u(i)+u(i-1));
    end
      figure(1);
      pause(0.5);
      hold on;
      plot(x,u);
end
 
