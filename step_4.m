%step 4
%1D Burgers equations 

clear all;
close all;
clc;

%constant
mu=0.3;
vis=mu;
pi=3.1416;
%Mesh 
xl=2*pi;
dx=0.1;
x=0:dx:xl;
size_x=size(x);
N=size_x(2);
%Initial condition
% phi=exp(-x.^2/(4*mu))+exp(-(x-2*pi).^2/(4*mu));
% dphi=((-0.5.*x./mu).*exp(-0.25.*x.^2)./mu)+...
%     (0.5.*((2*pi)-x)./mu).*exp(-0.25.*(((2*pi)-x.^2)./mu));

for i=1:N
     phi(i)=exp(-0.25*(x(i)^2)/vis)+exp(-0.25*(((2*pi)-x(i))^2)/vis);
     dphi(i)=(-0.5*x(i)/vis)*exp(-0.25*(x(i)^2)/vis)+(0.5*((2*pi)-x(i))/vis)*exp(-0.25*(((2*pi)-x(i))^2)/vis);
     u(i)=(-2*vis*(dphi(i)/phi(i)))+4;
end
% u=-(2*mu./phi).*(dphi)+4;




%time 
dt=0.01;
tf=1;
nt=round(tf/dt);

%boundary condition

u(1)=u(N); %periodic boundary conditions

%analytical solution


%main body

for t=1:nt
    u(1)=u(N);
    for i=2:N-1
        u(i)=u(i)-(dt/dx)*(u(i)*(u(i)-u(i-1)))+...
            ((mu*dt)/dx.^2)*(u(i+1)-2*u(i)+u(i-1));
    end
    figure (1);
    axis([0 6 0 6])
    pause(0.1);
    hold on;
    plot(x,u);
end
