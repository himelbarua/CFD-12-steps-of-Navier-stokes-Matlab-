%Loren barbaba
%step:1, 1D Linear convectiion
%-----------------------------------
clear all;
close all;
clc;

%Mesh 
dx=0.1;
xl=2;
N=round(xl/dx);
n1=round(0.5/dx);
n2=round(1/dx);
x=0:dx:xl;
%Initial value
size_x=size(x);
u=ones(size_x(2),1);
dt=0.01; %time step size
tf=1;%final time
c=1;
nt=round(tf/dt);
u(n1:n2)=2;

for t=1:nt
    for i=2:N
        %u(i)=u(i)-(c*dt/dx)*(u(u+1)-u(i));
        u(i)=u(i)-(c*dt/dx)*(u(i)-u(i-1));
    end
    figure(1);
    pause (0.5);
   plot(x,u);
   hold on;
end     
