%step 3
%2D Laplace equations
% Note: (mu*dt)/dx^2 has to be <1


clear all;
close all;
clc;

%constant
mu=0.3;

%Mesh
xl=2;yl=2;
x1=0.5;x2=1;
x2=1;
dx=0.025;dy=0.025;
Nx=round(xl/dx);Ny=round(yl/dy);
x=0:dx:xl;
y=0:dy:yl;

%Initial condition
dt=0.01;
tf=1;
nt=round(tf/dt);
size_x=size(x);size_y=size(y);
Nx=size_x(2);Ny=size_y(2);
p=zeros(Nx,Ny);
%Boundary condiitons

% p(1,1:Ny)=0;
% p(Nx,1:Ny)=Ny:-1:1;
% p(1:Nx,1)=p(1:Nx,2);
% p(1:Nx,Ny)=p(1:Nx,Ny-1);

  p(Nx,:)=y;
  p(2:Nx-1,1)=p(2:Nx-1,2);
  p(2:Nx-1,Ny)=p(2:Nx-1,Ny-1);
  pfinal=p;


%loop

for t=1:nt
    p=pfinal;
    for i=2:Nx-1
        for j=2:Ny-1
            pfinal(i,j)=(1/(dx.^2+dy.^2))*((dy.^2)*(p(i+1,j)+p(i-1,j))+(dx.^2)*(p(i,j+1)+p(i,j-1)));       
%             p(i,j) = ((pd(i+1,j)+pd(i-1,j))*dy^2+ (pd(i,j-1)+pd(i,j+1))*dx^2 )/(dx^2+dy^2)/2;
        end
    end
      figure(1);
      title('2D Laplace equation');
      pause(0.001);
      surf(x,y,pfinal);
       shading interp;
end
 
