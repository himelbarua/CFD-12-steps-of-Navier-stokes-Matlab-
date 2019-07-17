%2D poisson equations
%step_10

clear all;
close all;
clc;

%Mesh
xl=2;yl=2;
dx=0.1;dy=0.1;
x=0:dx:xl;y=0:dy:yl;
size_x=size(x);size_y=size(y);
Nx=size_x(2);Ny=size_y(2);

%Initialization
p=zeros(Nx,Ny);
pd=zeros(Nx,Ny);
b=zeros(Nx,Ny);
n_ite=5000;

%constant
b((Nx-1)*0.25,(Ny-1)*0.25)=100;
b((Nx-1)*0.75,(Ny-1)*0.75)=-100;

%Boundary conditions
p(1,:)=0;
p(Nx,:)=0;
p(:,1)=0;
p(:,Ny)=0;

for t=1:n_ite
    pd=p;
    for i=2:Nx-1
        for j=2:Ny-1
            p(i,j)=(((pd(i+1,j)+pd(i-1,j))*dy^2+(pd(i,j+1)+pd(i,j-1))*dx^2)...
                -b(i,j)*dx^2*dy^2)*(1/(2*(dx^2+dy^2)));
        end
    end
    
    
    figure(1);
    title('2D poisson equations');
    pause(0.01);
    surf(x,y,p);
    xlabel('x');ylabel('y');
end
            
            
            
            
            
            