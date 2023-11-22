lambda=13e-9;%wavelength
k=2*pi/lambda;%wavenumber
L=1e-4;%distance of image plane to source plane
dslitwidth=1e-7;% width between fringes
dwidth=8e-6;%width between slit
N=700;%number of slits
Ndx=1000;%number of point
E3a=2*1e-2;%amplitude of E3
x1=linspace(-dwidth/2,-dwidth/2-N*dslitwidth,N);%dimension of slit 1 at coordinate of source plane
x2=linspace(dwidth/2,dwidth/2+N*dslitwidth,N);%dimension of slit 2 at coordinate of source plane
x3=zeros(1,N);%dimension of position of external E field

dx1=L*lambda./(N*dwidth);% width
dx=linspace(-N*dx1,N*dx1,Ndx);% display width
phi1=zeros(Ndx,N);% phase  1
phi2=zeros(Ndx,N);% phase  2
phi3=zeros(Ndx,N);% phase  3
deltaphi=zeros(Ndx,N);% phase difference between 1 and 2
deltaphi2=zeros(Ndx,N);% phase difference between 2 and 3
deltaphi3=zeros(Ndx,N);% phase difference between 1 and 3
deltaphi4=zeros(Ndx,N);% phase difference between 1,2 and 3
% E1=zeros(Ndx,1);
% E2=zeros(Ndx,1);
E1=zeros(Ndx,N);% electric field 1
E2=zeros(Ndx,N);% electric field 2
E3=zeros(Ndx,N);% electric field 3 weak field
I1=zeros(Ndx,N);% intensity 1
I2=zeros(Ndx,N);% intensity 2
I3=zeros(Ndx,N);% intensity 3 weak field
for j=1:Ndx
for i=1:N
phi1(j,i)=2*pi*sqrt(L^2+(dx(j)-x1(i)).^2)./lambda;
phi2(j,i)=2*pi*sqrt(L^2+(dx(j)-x2(i)).^2)./lambda;
phi3(j,i)=2*pi*sqrt(L^2+(dx(j)-x3(i)).^2)./lambda;
deltaphi(j,i)=phi1(j,i)-phi2(j,i);
deltaphi2(j,i)=phi2(j,i)-phi3(j,i);
deltaphi3(j,i)=phi1(j,i)-phi3(j,i);
deltaphi4(j,i)=phi1(j,i)-phi2(j,i)-phi3(j,i);
E1(j,i)=(cos(k*(L)+(phi1(j,i))));
E2(j,i)=(cos(k*(L)+(phi2(j,i))));
E3(j,i)=E3a*(cos(k*(L)+(phi3(j,i))));
I1(j,i)=E1(j,i).*conj(E1(j,i));
I2(j,i)=E2(j,i).*conj(E2(j,i));
I3(j,i)=E3(j,i).*conj(E3(j,i));
end
end
I=zeros(Ndx,1);%intensity
Icorr=zeros(Ndx,1);%intensity with external E field

deltaphi1=zeros(Ndx,1);% total phase difference between 1 and 2
deltaphi22=zeros(Ndx,1);% total phase difference between 2 and 3
deltaphi33=zeros(Ndx,1);% total phase difference between 1 and 3
deltaphi44=zeros(Ndx,1);% total phase difference between 1,2 and 3
for i=1:Ndx
deltaphi1(i)=sum(deltaphi(i,:));
deltaphi22(i)=sum(deltaphi2(i,:));
deltaphi33(i)=sum(deltaphi3(i,:));
deltaphi44(i)=sum(deltaphi4(i,:));

I(i)=sum(I1(i,:)+I2(i,:)+2*sqrt(I1(i,:).*I2(i,:)).*cos(deltaphi(i,:)));
Icorr(i)=sum(I1(i,:)+I2(i,:)+2*sqrt(I1(i,:).*I2(i,:)).*cos(deltaphi(i,:))+2*sqrt(I2(i,:).*I3(i,:)).*cos(deltaphi22(i,:))+2*sqrt(I1(i,:).*I3(i,:)).*cos(deltaphi33(i,:))+sqrt(I1(i,:).*I2(i,:).*I3(i,:)).*cos(deltaphi44(i,:)));
% E1(i)=sum(cos(N*k*(L)+(phi1(i,:))));
% E2(i)=sum(cos(N*k*(L)+(phi2(i,:))));
end
E=E1+E2;
% I=E1.*conj(E1)+E2.*conj(E2)+2*sqrt(E1.*conj(E1)+E2.*conj(E2)).*cos(deltaphi1);
figure(1);
plot(dx,I,'r-');
figure(2);
plot(dx,Icorr,'r-');