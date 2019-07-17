%Loren barbaba
%step:1, 2D Linear convectiion
%-----------------------------------
clear all;
close all;
clc;

%Mesh 
dx=0.025;dy=dx;
xl=2;yl=2;
Nx=round(xl/dx);
Ny=round(yl/dy);
nx1=round(0.5/dx);
nx2=round(1/dx);
ny1=round(0.5/dy);
ny2=round(1/dy);
x=0:dx:xl;
y=0:dy:yl;
%Initial value
size_x=size(x);
size_y=size(y);
% u=ones(size_x(2),size_y(2));
u=zeros(size_x(2),size_y(2));
dt=0.001; %time step size
tf=0.5;%final time
c=1;
nt=round(tf/dt);
% for i=nx1:nx2
%     for j=ny1:ny2
%         u(i,j)=2;
%     end
% end
u(nx1:nx2,ny1:ny2)=2;

for t=1:nt
       for i=2:Nx
        for j=2:Ny
        %u(i)=u(i)-(c*dt/dx)*(u(u+1)-u(i));
        u(i,j)=u(i,j)-(c*dt/dx)*(u(i,j)-u(i-1,j))-(c*dt/dy)*(u(i,j)-u(i,j-1));
        end
    end
    figure(1);
    title('2D linear advection');
    axis([0 2 0 2 0 2]);
    F(t)=getframe(gcf);
    pause (0.001);
    surf(x,y,u);
    shading interp;
    
%    hold on;
end     

video=VideoWriter('2D linear advection', 'Uncompressed AVI');
open(video)
writeVideo(video,F)
close(video);