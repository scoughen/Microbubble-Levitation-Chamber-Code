clear
close all
clc

file = "LiLens500kHzHighResPhaseScan10SampleAvg0SecDelay3mmOffsetFromTransducerNewAmp9.6Vpp.csv";       
xRes = 0.5;
yRes = 1;
zRes = 0.5;
vToMPa = 0.8; %500kHz = 0.8V/MPa       2.25MHz = 0.92V/MPa

Layer = 5;

M = readmatrix(file);

x = M([1:end],1);
y = M([1:end],2);
z = M([1:end],3);
aV = M([1:end],4); %(V)
a = aV/vToMPa; %(MPa)
a = [a;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];
pha = (M(1:end,5) - (max(M(1:end,5))-180)) * (pi/180);
pha = [pha;30;30;30;30;30;30;30;30;30];

x1 = min(x):xRes:max(x);
% x1 = x1(1:end-1);
y1 = min(y):yRes:max(y);
z1 = min(z):zRes:max(z);

centerX = 94;
centerY = 165;
topZ = max(z);

x1 = (x1-centerX)/1000;
y1 = (y1-centerY)/1000;
z1 = -(z1-topZ-3)/1000;

xMax = find(x1==-min(x1));
xCrop = x1(1:xMax);

[X,Y,Z] = meshgrid(xCrop,y1,z1);

ptsPerLayer = length(x1)*length(y1);
numLayers = ceil(length(z)/ptsPerLayer);

startIndex = 1;

for layer = 1:numLayers
    for i = 1:length(y1)
        endIndex = startIndex+length(x1)-1;
        A(i,:,layer) = a(startIndex:endIndex);
        startIndex = endIndex+1;
        Pha(i,:,layer) = pha((i-1)*length(x1)+1 + (ptsPerLayer*(layer-1)):i*length(x1) + (ptsPerLayer*(layer-1)));
%         Pha(i,:,layer) = pha(startIndex:endIndex);
        if mod(i,2) == 0
            A(i,:,layer) = flip(A(i,:,layer));
            Pha(i,:,layer) = flip(Pha(i,:,layer));
        end  
    end
    A1 = imgaussfilt(A,1);
    A2 = imgaussfilt(A,2);
    
end

A = flip(A,3);
A1 = flip(A1,3);
A2 = flip(A2,3);
Pha = flip(Pha,3);


ACrop = A(:,1:length(xCrop),:);
A1Crop = A1(:,1:length(xCrop),:);
A2Crop = A2(:,1:length(xCrop),:);
PhaCrop = Pha(:,1:length(xCrop),:);

yMidIndex = find(y1==0);

xx = reshape(X(yMidIndex,:,:),length(xCrop),[]);
zz = reshape(Z(yMidIndex,:,:),length(xCrop),[]);
AA = reshape(ACrop(yMidIndex,:,:),length(xCrop),[]);
AA1 = reshape(A1Crop(yMidIndex,:,:),length(xCrop),[]);
AA2 = reshape(A2Crop(yMidIndex,:,:),length(xCrop),[]);
Pha = reshape(PhaCrop(yMidIndex,:,:),length(xCrop),[]);

TopAvgScan = sum(AA(:,9), 'all') / length(A(:,1))
AreaAvgScan = sum(AA, 'all') / (length(A(:,1))*length(A(1,:)))

figure
surf(xx,zz,Pha, 'edgecolor','none')
set(gca, 'YDir', 'reverse')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
colormap jet
colorbar
caxis([-pi pi])
colorbar('Ticks',[-pi,-pi/2,0,pi/2,pi],'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
% ylim([0.003, 0.053])
view(2)


figure
H = surf(xx,zz,AA);
set(gca, 'YDir', 'reverse')
axis equal
xlabel('X (mm)')
ylabel('Y (mm)')
zlabel('Z (mm)')
colorbar
% ylim([0.003, 0.053])
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





