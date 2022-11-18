%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script plots a given scan in a specified vertical plane. This
%  vertical plane is orthogonal to the primary vertical scan plotting
%  script
%
%  The parameters that need to be set are:
%    Scan parameter:
%      file = file name of the scan
%      xRes, yRes, zRes = the x, y, and z resolution of the scan
%      vToMPa = the sensitivity of the needle hydrophone used for scanning
%      centerX, centerY = are the X and Y coordinates of the center of the
%                         acoustic field
% 
%  S. Coughenour - Nov. 17, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% parameters
file = "LiLens500kHzHighResPhaseScanOrthogonal5SampleAvg0SecDelay3mmOffsetFromTransducerNewAmp9.6VppPart3.csv";       
xRes = 1;
yRes = 0.5;
zRes = 0.5;
vToMPa = 0.8; %500kHz = 0.8V/MPa       2.25MHz = 0.92V/MPa


%reading and processing scan data
M = readmatrix(file);

x = M(1:end,1);
y = M(1:end,2);
z = M(1:end,3);
aV = M(1:end,4); %(V)
a = aV/vToMPa; %(MPa)
% a = [a;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];
pha = (M(1:end,5) - (max(M(1:end,5))-180)) * (pi/180);
% pha = [pha;30;30;30;30;30;30;30;30;30];

x1 = min(x):xRes:max(x);
% x1 = x1(1:end-1);
y1 = min(y):yRes:max(y);
z1 = min(z):zRes:max(z);

centerX = 92;
centerY = 108;%105;
topZ = max(z);

x1 = (x1-centerX)/1000;
y1 = (y1-centerY)/1000;
z1 = -(z1-topZ-3)/1000;

yMax = find(y1==max(y1));%-min(y1));
yCrop = y1(1:yMax);

[X,Y,Z] = meshgrid(x1,yCrop,z1);

% rearrange scan data from vector into 2D matrix
ptsPerLayer = length(x1)*length(y1);
numLayers = ceil(length(z)/ptsPerLayer);
startIndex = 1;
A = zeros(length(x1), length(y1), length(z1));
Pha = zeros(length(x1), length(y1), length(z1));
for i = 1:length(a)
    A(x(i)*2-min(x*2)+1, y(i)*2-min(y*2)+1, z(i)*2-min(z*2)+1) = a(i);
    Pha(x(i)*2-min(x*2)+1, y(i)*2-min(y*2)+1, z(i)*2-min(z*2)+1) = pha(i);
end

for i = 2:11 %2:8 for part 1 and 2:6 for part 2 and 2:11 for part 3
    A(i,:,:) = [];
    Pha(i,:,:) = [];
end

A1 = imgaussfilt(A,1); %smoothing data
A2 = imgaussfilt(A,2);



% data manipulation of 2D matrix for plotting

% A = flip(A,3);
% A1 = flip(A1,3);
% A2 = flip(A2,3);
% Pha = flip(Pha,3);

ACrop = A(:,1:length(yCrop),:);
A1Crop = A1(:,1:length(yCrop),:);
A2Crop = A2(:,1:length(yCrop),:);
PhaCrop = Pha(:,1:length(yCrop),:);

xMidIndex = find(x1==0);

xx = reshape(Y(:,xMidIndex,:),length(yCrop),[]);
zz = reshape(Z(:,xMidIndex,:),length(yCrop),[]);
AA = reshape(ACrop(xMidIndex,:,:),length(yCrop),[]);
AA1 = reshape(A1Crop(xMidIndex,:,:),length(yCrop),[]);
AA2 = reshape(A2Crop(xMidIndex,:,:),length(yCrop),[]);
PPha = reshape(PhaCrop(xMidIndex,:,:),length(yCrop),[]);

TopAvgScan = sum(AA(:,9), 'all') / length(A(:,1))
AreaAvgScan = sum(AA, 'all') / (length(A(:,1))*length(A(1,:)))


%%% plot data

figure
surf(xx,zz,PPha, 'edgecolor','none')
set(gca, 'YDir', 'reverse')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
colormap jet
colorbar
caxis([-pi pi])
colorbar('Ticks',[-pi,-pi/2,0,pi/2,pi],'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
view(2)


figure
H = surf(xx,zz,AA);
set(gca, 'YDir', 'reverse')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
colorbar
set(H,'edgecolor','none')
view(0,90)

% figure
% H = surf(xx,zz,AA1);
% axis equal
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% colorbar
% set(H,'edgecolor','none')
% view(0,90)
% 
% figure
% H = surf(xx,zz,AA2);
% axis equal
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% colorbar
% set(H,'edgecolor','none')
% view(0,90)

% figure
% xslice = [ceil(median(x1))];
% yslice = [];
% zslice = [];
% s = slice(X,Y,Z,A1,xslice,yslice,zslice);
% xlabel('X (mm)')
% ylabel('Y (mm)')
% zlabel('Z (mm)')
% colorbar
% set(s,'edgecolor','none')
% view(90,0)





