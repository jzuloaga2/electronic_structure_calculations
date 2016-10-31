% This program plots the equilibrium density. 
clear;
r=12;cent=16;ang=0:0.01:2*pi;xp=r*cos(ang);yp=r*sin(ang);
myPara = load('myPara.txt','-ascii');
Nx=myPara(2);
Nz=myPara(3);
Lx=myPara(12);
Lz=myPara(13);
Wx=myPara(4);
Wz=myPara(5);

den = load('den.txt','-ascii');
pot= load('pot.txt','-ascii');

% Here we make the first (the original) part of the matrix:
ii=0;
denM = zeros(2*Nz+1,Nx);
potM = zeros(2*Nz+1,Nx);
for ix=1:Nx
    for iz = -Nz:Nz
      ii=ii+1;
      denM(iz+Nz+1,ix)=den(ii);
      potM(iz+Nz+1,ix)=pot(ii);
    end
end
% Now we expand the matrix to include the other half:
denM2 = zeros(2*Nz+1,2*Nx+1);
potM2 = zeros(2*Nz+1,2*Nx+1);
x = zeros(2*Nz+1,2*Nx+1);
z = zeros(2*Nz+1,2*Nx+1);
for ix=-Nx:Nx
    for iz = -Nz:Nz
    if(ix<0)
       denM2(iz+Nz+1,ix+Nx+1)=denM(iz+Nz+1,-ix);
       potM2(iz+Nz+1,ix+Nx+1)=potM(iz+Nz+1,-ix);
    elseif(ix==0)
       denM2(iz+Nz+1,ix+Nx+1)=denM(iz+Nz+1,1);
       potM2(iz+Nz+1,ix+Nx+1)=potM(iz+Nz+1,1);
    else
       denM2(iz+Nz+1,ix+Nx+1)=denM(iz+Nz+1,ix);
       potM2(iz+Nz+1,ix+Nx+1)=potM(iz+Nz+1,ix);
    end
      x(iz+Nz+1,ix+Nx+1)=ix*Wx;
      z(iz+Nz+1,ix+Nx+1)=iz*Wz;
    end
end
%Now we are ready to plot:
%plot outlines of MNP's


figure
mesh(x,z,denM2)
title('Density')
shading interp;
%axis equal;
%axis tight;
hold on;

figure
mesh(x,z,potM2)
title('Potential')
shading interp;
%axis equal;
%axis tight;
hold on;
