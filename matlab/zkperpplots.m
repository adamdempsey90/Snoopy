% This file reads in snoopy data and computes the FFT of the data
% so that we can view the k-spectrum - want z,kperp plots
% 
% Adrian Barker 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; %close all;
rep='../data/';
%rep='/projects/b1002/adrian/RB/128om10F1Lh1nu2p5kap2p5/';
%rep='/projects/b1002/adrian/RB/256om10F1Lh1nu4kap4/';
%rep='/projects/b1002/adrian/RB/HC128om30F1Lh0p3Lz1p4nu3p5kap3p5/';
%rep='/projects/b1002/adrian/RB/HC256om30F1Lh0p3Lz1p4nu4kap4/';
%rep='/projects/b1002/adrian/RB/HC128om30F1Lh0p6Lz1p4nu3p5kap3p5/';
%rep='/projects/b1002/adrian/RB/HC256om3F1Lh2Lz1p4nu3p5kap3p5/';
%rep='/projects/b1002/adrian/RB/HC128om20F1Lh0p7Lz1p4nu3p3kap3p3/';
%rep='/projects/b1002/adrian/RB/HC128om10th80F1Lh1Lz1p4nu3kap3/';

% frame number
number=0010; numfiles=1; numcontours = 50; khzplots=0;

for tindx=1:numfiles;
%Read VTK files
[V,time]=readVTK([rep,'v',num2str(number,'%0.4d'),'.vtk']);
t(tindx) = time;

% Fourier components
n = (size(V.vx,1)); vkx = fftn(V.vx)./n^3;
vkxz = vkx(1,1,:);
n = (size(V.vy,1)); vky = fftn(V.vy)./n^3;
vkyz = vky(1,1,:);
n = (size(V.vz,1)); vkz = fftn(V.vz)./n^3;
vkzz = vkz(1,1,:);
K(:,:,:,tindx) = 0.5*(vkx.*conj(vkx) + vky.*conj(vky) + vkz.*conj(vkz));
n = (size(V.th,1)); thk = fftn(V.th)./n^3;
P(:,:,:,tindx) = 0.5*(thk.*conj(thk));
n = (size(V.th,1)); uzthk = fftn(V.th.*V.vz)./n^3;
uzth(:,:,:,tindx) = uzthk(1,1,:);

%z,kperp spectra
n = size(V.vx,1);
for i=1:size(V.vz,3); %perform fft first, then abs them to avoid convolutions
vxkz(:,:,i,tindx)=fft2(V.vx(:,:,i))./n^2;
vykz(:,:,i,tindx)=fft2(V.vy(:,:,i))./n^2;
vzkz(:,:,i,tindx)=fft2(V.vz(:,:,i))./n^2;
thkz(:,:,i,tindx)=fft2(V.th(:,:,i))./n^2;
uz2spec(:,:,i,tindx)=0.5*(vzkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx)));
uh2spec(:,:,i,tindx)=0.5*(vxkz(:,:,i,tindx).*conj(vxkz(:,:,i,tindx))+vykz(:,:,i,tindx).*conj(vykz(:,:,i,tindx)));
Kspec(:,:,i,tindx)=0.5*(vxkz(:,:,i,tindx).*conj(vxkz(:,:,i,tindx)) + vykz(:,:,i,tindx).*conj(vykz(:,:,i,tindx)) + vzkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx)));
Pspec(:,:,i,tindx)=0.5*thkz(:,:,i,tindx).*conj(thkz(:,:,i,tindx));
uzthspec(:,:,i,tindx)=conj(thkz(:,:,i,tindx)).*vzkz(:,:,i,tindx)+thkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx));
end;
number=number+1;
end;

%z,kperp plots
%now create z,kh plots
Kavspec = mean(abs(Kspec),4); nal=round(n/3);
Kavspeckh = zeros(size(Kspec,3),nal);
uz2avspec = mean(abs(uz2spec),4);
uz2avspeckh = zeros(size(uz2spec,3),nal);
uh2avspec = mean(abs(uh2spec),4);
uh2avspeckh = zeros(size(uh2spec,3),nal);
Pavspec = mean(abs(Pspec),4);
Pavspeckh = zeros(size(Pspec,3),nal);
uzthavspec = mean(abs(uzthspec),4);
uzthavspeckh = zeros(size(uzthspec,3),nal);
for m=1:nal;
number = 0;
 for i=1:size(Kspec,1);
  for j=1:size(Kspec,2);
   if(i^2+j^2 < (m+0.5)^2 && i^2+j^2 >= (m-0.5)^2)
    Kavspeckh(:,m) = Kavspeckh(:,m) + squeeze(Kavspec(i,j,:));
    uz2avspeckh(:,m) = uz2avspeckh(:,m) + squeeze(uz2avspec(i,j,:));
    uh2avspeckh(:,m) = uh2avspeckh(:,m) + squeeze(uh2avspec(i,j,:));
    Pavspeckh(:,m) = Pavspeckh(:,m) + squeeze(Pavspec(i,j,:));
    uzthavspeckh(:,m) = uzthavspeckh(:,m) + squeeze(uzthavspec(i,j,:));
    number = number + 1;
   end;
  end;
 end;
Kavspeckh(:,m) = Kavspeckh(:,m)./number; %mean
uz2avspeckh(:,m) = uz2avspeckh(:,m)./number; %mean
uh2avspeckh(:,m) = uh2avspeckh(:,m)./number; %mean
uzthavspeckh(:,m) = uzthavspeckh(:,m)./number; %mean
Pavspeckh(:,m) = Pavspeckh(:,m)./number; %mean
end;
z=(0:size(Kspec,3)-1)./size(Kspec,3);
kh=0:nal-1;
khmat=repmat(2*pi*kh.^2,size(z,2),1);

%1D horizontal spectra
%indexbl=10;
K1D=sum(Kavspeckh(:,:),1);
uz21D=sum(uz2avspeckh(:,:),1);
uh21D=sum(uh2avspeckh(:,:),1);
P1D=sum(Pavspeckh(:,:),1);
uzth1D=sum(uzthavspeckh(:,:),1);
figure(1);
loglog(kh,K1D.*kh.^2,kh,uz21D.*kh.^2,kh,uh21D.*kh.^2,kh,P1D.*kh.^2,kh,uzth1D.*kh.^2); hold on;
%loglog(kh,K1D.*kh.^2,'--',kh,uz21D.*kh.^2,'--',kh,uh21D.*kh.^2,'--',kh,P1D.*kh.^2,'--',kh,uzth1D.*kh.^2,'--'); hold on;
%loglog(kh,K1D.*kh.^2,'-.',kh,uz21D.*kh.^2,'-.',kh,uh21D.*kh.^2,'-.',kh,P1D.*kh.^2,'-.',kh,uzth1D.*kh.^2,'-.'); hold on;
loglog(kh,1e3*kh.^(-5/3),'k--',kh,1e3*kh.^(-2/3),'--');
y1=-6:0.01:3; x1=nal-1; loglog(x1,10.^(y1),'k')
ylim([1e-4,1e3]);
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
legend('K','uz^2','uh^2','b^2','uzb');

Ekperp=K1D.*kh; uzbkperp=uzth1D.*kh; Bkperp=P1D.*kh;
disp('kpeak |uk|^2:'); trapz(kh,Ekperp.*kh)./trapz(kh,Ekperp)
disp('kpeak |uzb_k|:'); trapz(kh,uzbkperp.*kh)./trapz(kh,uzbkperp)
disp('kpeak |b_k|^2:'); trapz(kh,Bkperp.*kh)./trapz(kh,Bkperp)
disp('kpeak |u_z,k|^2:'); trapz(kh,uz21D.*kh.^2)./trapz(kh,uz21D.*kh)








if(khzplots==1)
figure(2);
contour(kh,z,Kavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Kinetic spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

figure(3);
contour(kh,z,uz2avspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Vertical kinetic energy spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

figure(4);
contour(kh,z,uh2avspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Horizontal kinetic energy spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

figure(5);
contour(kh,z,Pavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Thermal spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

figure(6);
contour(kh,z,uzthavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Flux spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);
end
