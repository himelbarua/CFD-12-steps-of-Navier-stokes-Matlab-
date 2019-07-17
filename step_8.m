%step 4
%2D Burgers equations 

clear all;
close all;
clc;

%constant
mu=0.001;
pi=3.1416;

%Mesh 

xl=2;yl=2;
dx=0.05;dy=0.05;
x=0:dx:xl;y=0:dy:yl;
size_x=size(x);size_y=size(y);
Nx=size_x(2);Ny=size_y(2);
x1=0.5,x2=1;y1=0.5;y2=1;
nx1=round(x1/dx);
nx2=round(x2/dx);
ny1=round(y1/dy);
ny2=round(y2/dy);

%Initialization

u=ones(Nx,Ny);
v=ones(Nx,Ny);
u(nx1:nx2,ny1:ny2)=2;
v(nx1:nx2,ny1:ny2)=2;

%time 
dt=0.001;
tf=1;
nt=round(tf/dt);


%main body

for t=1:nt
    
    for i=2:Nx-1
        
        for j=2:Ny-1
            
        u(i,j)=u(i,j)-(dt/dx)*(u(i,j)*(u(i,j)-u(i-1,j)))...
            -(dt/dy)*v(i,j)*(u(i,j)-u(i,j-1))...
            +((mu*dt)/dx.^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))...
            +((mu*dt)/dy.^2)*(u(i,j+1)-2*u(i,j)+u(i,j-1));
        
        v(i,j)=v(i,j)-(dt/dx)*(u(i,j)*(v(i,j)-v(i-1,j)))...
            -(dt/dy)*(v(i,j)*(v(i,j)-v(i,j-1)))...
            +((mu*dt)/dx.^2)*(v(i+1,j)-2*v(i,j)+v(i-1,j))...
            +((mu*dt)/dy.^2)*(v(i,j+1)-2*(v(i,j)+v(i,j-1)));
        end
    end
    
    figure (2);
    shading interp;
    axis([0 2 0 2 0 2 ]);
    pause(0.001);
   surf(x,y,u);
    
    figure(3);
    shading interp;
    axis([0 2 0 2 0 2]);
    pause(0.001);
    surf(x,y,v);
    
end
