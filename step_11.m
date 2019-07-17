%lid driven cavity flow for incompressible flow
%Step 11

clear all;
close all;
clc;


%Mesh
xl=1;yl=1;
dx=0.05;dy=0.05;
x=0:dx:xl;y=0:dy:yl;
size_x=size(x);size_y=size(y);
Nx=size_x(2);Ny=size_y(2);

%constnat
mu=0.3;
rho=0.1;

%Initialization

u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
p=zeros(Nx,Ny);
pd=zeros(Nx,Ny);
b=zeros(Nx,Ny);

%time
tf=0.1;
dt=0.001;
Nt=round(tf/dt);
n_iter=5000;


 %Boundary conditions

 u(:,1)=0;
 u(1:Nx-1,Ny)=1;
 u(1,1:Ny-1)=0;
 u(Nx,1:Ny-1)=0;
 v(1,:)=0;
 v(Nx,:)=0;
 v(:,1)=0;
 v(:,Ny)=0;


%corner boundary conditions
p(1,1)=(p(2,1)+p(1,2))*0.5;
p(Nx,1)=(p(Nx-1,1)+p(Nx,2))*0.5;
p(1,Ny)=(p(1,Ny-1)+p(2,Ny))*0.5;
p(Nx,Ny)=(p(Nx-1,Ny)+p(Nx,Ny-1))*0.5;

u(1,Ny)=0.5*(u(1,Ny-1)+u(2,Ny));
u(Nx,Ny)=0.5*(u(Nx-1,Ny)+u(Nx,Ny-1));

%main body

for t=1:Nt 
   
for i=2:Nx-1
   for j=2:Ny-1
                      
%   b(i,j)=((u(i+1,j)-u(i-1,j)).^2*(0.25)*(1/dx^2)-0.5*(u(i,j+1)-u(i,j-1))*(v(i+1,j)-v(i-1,j))*(1/(dx*dy))...
%  -(v(i,j+1)-v(i,j-1)).^2*(0.25)*(1/dy^2))*(1/dt);


 b(i,j)=rho*(1/dt*((u(i+1,j)-u(i-1,j))/(2*dx)+(v(i,j+1)-v(i,j-1))/(2*dy))-...
     ((u(i+1,j)-u(i-1,j))/(2*dx))^2-2*((u(i,j+1)-u(i,j-1))/(2*dy)*(v(i+1,j)...
     -v(i-1,j))/(2*dx))-((v(i,j+1)-v(i,j-1))/(2*dy))^2);
          
    end
end    
%poisson equation solving      
 for ite=1:n_iter
       
     pd=p;
       
       for i=2:Nx-1
           for j=2:Ny-1
               
               p(i,j)=(pd(i+1,j)+pd(i-1,j))*(dy^2)+(pd(i,j+1)+pd(i,j-1))*(dx)^2*0.5*(1/(dx^2+dy^2))...
                   -0.5*rho*dx^2*dy^2*(1/(dx^2+dy^2))*b(i,j);
           end
       end 
       
 p(1,2:Nx-1)=p(2,2:Nx-1); 
 p(Ny,2:Nx-1)=0;
 p(2:Ny-1,1)=p(2:Ny-1,2);
 p(2:Ny-1,Nx)=p(2:Ny-1,Nx-1); 
    
 end
            
 for i=2:Nx-1
        
    for j=2:Ny-1
            
       u(i,j)=u(i,j)-(dt/dx)*u(i,j)*(u(i,j)-u(i-1,j))-v(i,j)*(dt/dy)*(u(i,j)-u(i,j-1))...
                -(dt/(2*rho*dx))*(p(i+1,j)-p(i-1,j))...
                +mu*((dt/dx^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))...
                +(dt/dy^2)*(u(i,j+1)-2*u(i,j)+u(i,j-1)));
            
        v(i,j)=v(i,j)-(dt/dx)*u(i,j)*(v(i,j)-v(i-1,j))-v(i,j)*(dt/dy)*(v(i,j)-v(i,j-1))...
                -(dt/(2*rho*dy))*(p(i+1,j)-p(i-1,j))...
                +mu*((dt/dx^2)*(v(i+1,j)-2*v(i,j)+v(i-1,j))...
                +(dt/dy^2)*(v(i,j+1)-2*v(i,j)+v(i,j-1)));               
            
        end
        
 end

 
    figure(1);
    pause(0.001);
    contourf(x,y,u);
%     surf(x,y,u);
%     
%     figure(2);
%     pause(0.001);
%     surf(x,y,p);
    
end

% figure(3);
% contour(x,y,u);
% 
% figure(4);
% contour(x,y,v);

    
    
            
