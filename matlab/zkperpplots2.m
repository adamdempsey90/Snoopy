%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in snoopy VTK datafile(s) and computes the DFT of the data
% to produce the kh spectrum as function of z
% 
% This is not fast when NX,NZ>256...
%
% Adrian Barker 2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all; close all;
rep='../data/';

% frame number
number=0010; numfiles=1; numcontours = 50;

for tindx=1:numfiles;
%Read VTK files
[V,time]=readVTK([rep,'v',num2str(number,'%0.4d'),'.vtk']);
t(tindx) = time;

% Compute DFT
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
uz2spec(:,:,i,tindx)=real(vzkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx)));
uh2spec(:,:,i,tindx)=real(vxkz(:,:,i,tindx).*conj(vxkz(:,:,i,tindx))+vykz(:,:,i,tindx).*conj(vykz(:,:,i,tindx)));
Kspec(:,:,i,tindx)=real(vxkz(:,:,i,tindx).*conj(vxkz(:,:,i,tindx)) + vykz(:,:,i,tindx).*conj(vykz(:,:,i,tindx)) + vzkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx)));
Pspec(:,:,i,tindx)=real(thkz(:,:,i,tindx).*conj(thkz(:,:,i,tindx)));
uzthspec(:,:,i,tindx)=real(thkz(:,:,i,tindx).*conj(vzkz(:,:,i,tindx)));
end;
number=number+1;
end;

%z,kperp plots
%now create z,kh plots
Kavspec = mean(abs(Kspec),4); 

nal = ceil(sqrt((size(V.vx,1)/2-1)^2+(size(V.vy,1)/2-1)^2));

%nal=round(n/3); %aliasing limit
Kavspeckh = zeros(size(Kspec,3),nal);
for m=1:nal;
number = 0;
 for i=1:size(Kspec,1);
  for j=1:size(Kspec,2);
   if(i^2+j^2 < (m+0.5)^2 && i^2+j^2 >= (m-0.5)^2)
    Kavspeckh(:,m) = Kavspeckh(:,m) + squeeze(Kavspec(i,j,:));
    number = number + 1;
   end;
  end;
 end;
Kavspeckh(:,m) = Kavspeckh(:,m)./number; %mean
end;
figure
z=(0:size(Kspec,3)-1)./size(Kspec,3);
kh=0:nal-1;
khmat=repmat(2*pi*kh.^2,size(z,2),1);
contour(kh,z,Kavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Kinetic spectra cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

uz2avspec = mean(abs(uz2spec),4);
uz2avspeckh = zeros(size(uz2spec,3),nal);
for m=1:nal;
number = 0;
 for i=1:size(uz2spec,1);
  for j=1:size(uz2spec,2);
	local_kx = mod((i-1) + size(V.vx,1)/2,size(V.vx,1)) - size(V.vx,1)/2;
	local_ky = mod((j-1) + size(V.vx,2)/2,size(V.vx,2)) - size(V.vx,2)/2;
	local_kh = sqrt(local_kx^2 + local_ky^2);
	if abs(local_kh) < abs(m-1+.5) && abs(local_kh) >= abs(m-1-.5)

%   if((i-1)^2+(j-1)^2 < (m-1+0.5)^2 && (i-1)^2+(j-1)^2 >= (m-1-0.5)^2)
%    m = floor(sqrt((i-1)^2+(j-1)^2)+.5)
    uz2avspeckh(:,m) = uz2avspeckh(:,m) + squeeze(uz2avspec(i,j,:));
    number = number + 1;
   end;
  end;
 end;
%uz2avspeckh(:,m) = uz2avspeckh(:,m); %mean
end;
figure;
contour(kh,z,uz2avspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Vertical kinetic energy spectrum cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

uh2avspec = mean(abs(uh2spec),4);
uh2avspeckh = zeros(size(uh2spec,3),nal);
for m=1:nal;
number = 0;
 for i=1:size(uh2spec,1);
  for j=1:size(uh2spec,2);
   if(i^2+j^2 < (m+0.5)^2 && i^2+j^2 >= (m-0.5)^2)
    uh2avspeckh(:,m) = uh2avspeckh(:,m) + squeeze(uh2avspec(i,j,:));
    number = number + 1;
   end;
  end;
 end;
uh2avspeckh(:,m) = uh2avspeckh(:,m)./number; %mean
end;
figure;
contour(kh,z,uh2avspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Horizontal kinetic energy spectrum cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

Pavspec = mean(abs(Pspec),4);
Pavspeckh = zeros(size(Pspec,3),nal);
for m=1:nal;
number = 0;
   for i=1:size(Pspec,1);
     for j=1:size(Pspec,2);
        if(i^2+j^2 < (m+0.5)^2 && i^2+j^2 >= (m-0.5)^2)
	  Pavspeckh(:,m) = Pavspeckh(:,m) + squeeze(Pavspec(i,j,:));
	  number = number + 1;
	end;
     end;
   end;
Pavspeckh(:,m) = Pavspeckh(:,m)./number; %mean
end;
figure;
contour(kh,z,Pavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Thermal spectrum cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

uzthavspec = mean(abs(uzthspec),4);
uzthavspeckh = zeros(size(uzthspec,3),nal);
for m=1:nal;
number = 0;
   for i=1:size(uzthspec,1);
     for j=1:size(uzthspec,2);
	local_kx = mod((i-1) + size(V.vx,1)/2,size(V.vx,1)) - size(V.vx,1)/2;
	local_ky = mod((j-1) + size(V.vx,2)/2,size(V.vx,2)) - size(V.vx,2)/2;
	local_kh = sqrt(local_kx^2 + local_ky^2);
	if abs(local_kh) < abs(m-1+.5) && abs(local_kh) >= abs(m-1-.5)
   %     if((i-1)^2+(j-1)^2 < (m-1+0.5)^2 && (i-1)^2+(j-1)^2 >= (m-1-0.5)^2)
	  uzthavspeckh(:,m) = uzthavspeckh(:,m) + squeeze(uzthavspec(i,j,:));
	  number = number + 1;
	end;
     end;
   end;
%uzthavspeckh(:,m) = uzthavspeckh(:,m); %mean
end;
figure;
contour(kh,z,uzthavspeckh.*khmat,numcontours)
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
ylabel('z','fontsize',18);
title('Flux spectrum cylindrically averaged vs z','fontsize',18);
colorbar;
ylim([0 1]);
xlim([0 nal]);

%1D horizontal spectrum
%indexbl=20; %to remove heating/cooling and boundary layers
indexs=intersect(find(V.z >=-.25),find(V.z <= .25));
K1D=sum(Kavspeckh(indexs,:),1)./numel(indexs);
uz21D=sum(uz2avspeckh(indexs,:),1)./numel(indexs);
uh21D=sum(uh2avspeckh(indexs,:),1)./numel(indexs);
P1D=sum(Pavspeckh(indexs,:),1)./numel(indexs);
uzth1D=sum(uzthavspeckh(indexs,:),1)./numel(indexs);
figure; plot(kh,uz21D,'x'); title('matlab'); xlim([0 10]);
figure;
loglog(kh,K1D.*kh.^2,kh,uz21D.*kh.^2,kh,uh21D.*kh.^2,kh,P1D.*kh.^2,kh,uzth1D.*kh.^2); hold on;
loglog(kh,10.0*kh.^(-5/3),'k--',kh,10.0*kh.^(-2/3),'--');
set(gca,'fontsize',16)
xlabel('kh','fontsize',18);
legend('K','uz^2','uh^2','b^2','uzb');
