% Incompressible flow inside channel
%Step 11

clear all;
close all;
clc;

%Mesh
xl=2;yl=2;
dx=0.02;dy=0.02;
x=0:dx:xl;y=0:dy:yl;
size_x=size(x);size_y=size(y);
Nx=size_x(2);Ny=size_y(2);

%constnat
mu=0.3;
rho=0.1;
F=1;

%Initialization

u=zeros(Nx,Ny);
v=zeros(Nx,Ny);
p=zeros(Nx,Ny);
pd=zeros(Nx,Ny);
b=zeros(Nx,Ny);

%time
tf=0.1;
dt=0.0001;
Nt=round(tf/dt);
n_iter=5000;

%Boundary conditions

u(1:Nx,1)=0;
u(1:Nx,Ny)=0;
v(1:Nx,1)=0;
v(1:Nx,Ny)=0;

p(1:Nx,1)=p(1:Nx,2);
p(1:Nx,Ny)=p(1:Nx,Ny-1);


%corner boundary conditions
% p(1,1)=(p(2,1)+p(1,2))*0.5;
% p(Nx,1)=(p(Nx-1,1)+p(Nx,2))*0.5;
% p(1,Ny)=(p(1,Ny-1)+p(2,Ny))*0.5;
% p(Nx,Ny)=(p(Nx-1,Ny)+p(Nx,Ny-1))*0.5;
% 
% u(1,Ny)=0.5*(u(1,Ny-1)+u(2,Ny));
% u(Nx,Ny)=0.5*(u(Nx-1,Ny)+u(Nx,Ny-1));

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

%Boundary condition for b
%@ x=0
for j=2:Ny-1

b(1,j)=rho*(1/dt*((u(2,j)-u(Nx,j))/(2*dx)+(v(1,j+1)-v(1,j-1))/(2*dy))-...
     ((u(2,j)-u(Nx,j))/(2*dx))^2-2*((u(1,j+1)-u(1,j-1))/(2*dy)*(v(2,j)...
     -v(Nx,j))/(2*dx))-((v(1,j+1)-v(1,j-1))/(2*dy))^2);
end

 %@x=2
  for j=2:Ny-1
     
 b(1,j)=rho*(1/dt*((u(2,j)-u(Nx,j))/(2*dx)+(v(1,j+1)-v(1,j-1))/(2*dy))-...
     ((u(2,j)-u(Nx,j))/(2*dx))^2-2*((u(1,j+1)-u(1,j-1))/(2*dy)*(v(2,j)...
     -v(Nx,j))/(2*dx))-((v(1,j+1)-v(1,j-1))/(2*dy))^2);
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
       
  %Boundary condition for Pressure
  %@x=0
  
  for j=2:Ny-1 
  p(1,j)=(pd(2,j)+pd(Nx,j))*(dy^2)+(pd(1,j+1)+pd(1,j-1))*(dx)^2*0.5*(1/(dx^2+dy^2))...
                   -0.5*rho*dx^2*dy^2*(1/(dx^2+dy^2))*b(1,j);  
  end
  
  %@x=Nx
  for j=2:Ny-1
  
   p(Nx,j)=(pd(1,j)+pd(Nx-1,j))*(dy^2)+(pd(Nx,j+1)+pd(Nx,j-1))*(dx)^2*0.5*(1/(dx^2+dy^2))...
                   -0.5*rho*dx^2*dy^2*(1/(dx^2+dy^2))*b(Nx,j);
  end
  
      
p(1:Nx,1)=p(1:Nx,2);
p(1:Nx,Ny)=p(1:Nx,Ny-1);

 end
            
 for i=2:Nx-1
        
    for j=2:Ny-1
            
       u(i,j)=u(i,j)-(dt/dx)*u(i,j)*(u(i,j)-u(i-1,j))-v(i,j)*(dt/dy)*(u(i,j)-u(i,j-1))...
                -(dt/(2*rho*dx))*(p(i+1,j)-p(i-1,j))...
                +mu*((dt/dx^2)*(u(i+1,j)-2*u(i,j)+u(i-1,j))...
                +(dt/dy^2)*(u(i,j+1)-2*u(i,j)+u(i,j-1)))...
                +F*dt;
            
        v(i,j)=v(i,j)-(dt/dx)*u(i,j)*(v(i,j)-v(i-1,j))-v(i,j)*(dt/dy)*(v(i,j)-v(i,j-1))...
                -(dt/(2*rho*dy))*(p(i+1,j)-p(i-1,j))...
                +mu*((dt/dx^2)*(v(i+1,j)-2*v(i,j)+v(i-1,j))...
                +(dt/dy^2)*(v(i,j+1)-2*v(i,j)+v(i,j-1)));               
            
        end
        
 end
 
 %Boundary condition for u
 %@x=0
 for j=2:Ny-1
     
 u(1,j)=u(1,j)-(dt/dx)*u(1,j)*(u(1,j)-u(Nx,j))-v(1,j)*(dt/dy)*(u(1,j)-u(1,j-1))...
                -(dt/(2*rho*dx))*(p(2,j)-p(Nx,j))...
                +mu*((dt/dx^2)*(u(2,j)-2*u(1,j)+u(Nx,j))...
                +(dt/dy^2)*(u(1,j+1)-2*u(1,j)+u(1,j-1)))...
                +F*dt;
 end
 
 %@x=2
 
 for j=2:Ny-1
     
      u(Nx,j)=u(Nx,j)-(dt/dx)*u(Nx,j)*(u(Nx,j)-u(Nx-1,j))-v(Nx,j)*(dt/dy)*(u(Nx,j)-u(Nx,j-1))...
                -(dt/(2*rho*dx))*(p(1,j)-p(Nx-1,j))...
                +mu*((dt/dx^2)*(u(1,j)-2*u(Nx,j)+u(Nx-1,j))...
                +(dt/dy^2)*(u(Nx,j+1)-2*u(Nx,j)+u(Nx,j-1)))...
                +F*dt;
 end
 
     
 %Boundary condiiton for v
 %@x=0
 for j=2:Ny-1
     
     v(1,j)=v(1,j)-(dt/dx)*u(1,j)*(v(1,j)-v(Nx,j))-v(1,j)*(dt/dy)*(v(1,j)-v(1,j-1))...
                -(dt/(2*rho*dy))*(p(2,j)-p(Nx,j))...
                +mu*((dt/dx^2)*(v(2,j)-2*v(1,j)+v(Nx,j))...
                +(dt/dy^2)*(v(1,j+1)-2*v(1,j)+v(1,j-1))); 
 end
 
 %@x=2
 
 for j=2:Ny-1
     
 v(Nx,j)=v(Nx,j)-(dt/dx)*u(Nx,j)*(v(Nx,j)-v(Nx-1,j))-v(Nx,j)*(dt/dy)*(v(Nx,j)-v(Nx,j-1))...
                -(dt/(2*rho*dy))*(p(1,j)-p(Nx-1,j))...
                +mu*((dt/dx^2)*(v(1,j)-2*v(Nx,j)+v(Nx-1,j))...
                +(dt/dy^2)*(v(Nx,j+1)-2*v(Nx,j)+v(Nx,j-1)));
 end
 
 
 %Other Boundary conditions

u(1:Nx,1)=0;
u(1:Nx,Ny)=0;
v(1:Nx,1)=0;
v(1:Nx,Ny)=0;

p(1:Nx,1)=p(1:Nx,2);
p(1:Nx,Ny)=p(1:Nx,Ny-1);
 
 
 
%  figure(1);
%  pause(0.001);
%  contourf(x,y,u);
%  figure (2);
%  quiver(x,y,u,v);
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
 figure (2);
 quiver(x,y,u.',v.',2)    
%     
            
