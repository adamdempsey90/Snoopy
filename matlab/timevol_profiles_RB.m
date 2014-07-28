%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file reads in a SNOOPY profiles.dat file and computes
% various quantities for the Rotating Convection problem
%
% Adrian Barker 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
nz=128; %This should be NZ/2 from code = number of vertical grid points
kappa=10^(-4); %conductivity for conductive flux
rep='~/snoopyRBAJB/data/';
mat=load([rep,'profiles.dat'],'-ascii');
z1=reshape(mat(:,1),nz,size(mat,1)/nz); z=z1(:,1);
uxrms=reshape(sqrt(mat(:,2)),nz,size(mat,1)/nz);
uyrms=reshape(sqrt(mat(:,3)),nz,size(mat,1)/nz);
uhrms=reshape(sqrt(mat(:,2)+mat(:,3)),nz,size(mat,1)/nz);
urms=reshape(sqrt(mat(:,2)+mat(:,3)+mat(:,4)),nz,size(mat,1)/nz);
uzrms=reshape(sqrt(mat(:,4)),nz,size(mat,1)/nz);
th2=reshape(mat(:,5),nz,size(mat,1)/nz);
th2sqrt=sqrt(th2);
temp1=reshape(mat(:,6),nz,size(mat,1)/nz);
temp=temp1(:,1:2:end); %take out second half due to z=1->2,rather than 0->1
thvz=reshape(mat(:,7),nz,size(mat,1)/nz);
th=reshape(mat(:,8),nz,size(mat,1)/nz);
ux=reshape(mat(:,9),nz,size(mat,1)/nz);
uy=reshape(mat(:,10),nz,size(mat,1)/nz);
thrms=sqrt(th2-th.^2);
dzthtot=reshape(mat(:,11),nz,size(mat,1)/nz);%total including background
flux=thvz-kappa*dzthtot; %total flux
for i=1:size(temp,2);
dz=z(2)-z(1);
dbdz(:,i)=gradient(th(:,i),dz); %using Matlab
db2dz2(:,i)=gradient(dbdz(:,i),dz);
end;
tsteps=1:size(temp,2);
time=0.2*tsteps; %physical time

mean_start=int16(size(temp,2)/2); %start averages at 1/2 way point

figure(1);
plot(z,mean(temp(:,(mean_start/2):end),2),'x'); hold on;
set(gca,'fontsize',16);
ylabel('T','fontsize',18);
xlabel('z','fontsize',18);
indxbulk=int16(nz/3);
[polybulk,Sbulk] = polyfit(z(indxbulk:end-indxbulk),mean(temp(indxbulk:end-indxbulk,(mean_start/2):end),2),1);
fnpoly = polyval(polybulk,z);
plot(z,fnpoly,'r');
N2=polybulk(1);
N=sqrt(abs(N2));

figure(2);
for timeind=1:size(temp,2);
poly=polyfit(z(indxbulk:end-indxbulk),temp(indxbulk:end-indxbulk,timeind),1);
polys(timeind)=poly(1);
end;
plot(time,polys);
set(gca,'fontsize',16);
ylabel('N^2_bulk','fontsize',18);
xlabel('t','fontsize',18);
disp('Bulk slope:'); N2 %from mean T profile
disp('RMS error'); sqrt(mean((polys(mean_start:end)-mean(polys(mean_start:end))).^2))
disp('uz:'); mean(uzrms(end/2,mean_start:end))
disp('RMS uz:'); sqrt(mean((uzrms(end/2,mean_start:end)-mean(uzrms(end/2,mean_start:end))).^2))
disp('u:'); mean(urms(end/2,mean_start:end))
disp('th:'); mean(thrms(end/2,mean_start:end))
disp('F:'); mean(flux(end/2,mean_start:end))

figure(3);
plot(z,mean(thvz(:,mean_start:end),2),z,mean(thrms(:,mean_start:end).*uzrms(:,mean_start:end),2),z,mean(-kappa*dbdz(:,mean_start:end),2),z,mean(uzrms(:,mean_start:end).^2,2)*N,z,mean(thrms(:,mean_start:end).^2,2)./N,z,-kappa.*mean(dzthtot(:,mean_start:end),2))
set(gca,'fontsize',16);
ylabel('Flux','fontsize',18);
xlabel('z','fontsize',18);
legend('<thuz>','<thrms><uzrms>','-kapdbdz','uz^2N','deltab^2/N','kapdzthtot');

figure(4);
plot(z,mean(uzrms(:,mean_start:end),2),z,mean(uhrms(:,mean_start:end),2),z,mean(urms(:,mean_start:end),2))
set(gca,'fontsize',16);
ylabel('u','fontsize',18);
xlabel('z','fontsize',18);
legend('uz','uh','u');

figure(5);
plot(z,mean(thrms(:,mean_start:end),2),z,mean(uzrms(:,mean_start:end),2).*N,z,N.*z./z)
legend('deltab','deltab~uzN','deltab~ N^2 L')

figure(6);
plot(time,uzrms(end/2,1:2:end))
ylabel('uz','fontsize',18);
xlabel('t','fontsize',18);
