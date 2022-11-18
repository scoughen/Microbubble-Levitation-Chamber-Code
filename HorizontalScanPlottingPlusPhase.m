%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script plots a given scan of a series of horizontal planes
%
%  The parameters that need to be set are:
%    Scan parameter:
%      file = file name of the scan
%      xRes, yRes, zRes = the x, y, and z resolution of the scan
%      vToMPa = the sensitivity of the needle hydrophone used for scanning
% 
%  S. Coughenour - Nov. 17, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

% parameters
file = "LiLens500kHzHighResPhaseScanOrthogonal5SampleAvg0SecDelay30mmOffsetFromTransducerNewAmp9.6Vpp.csv";
xRes = 0.5;
yRes = 0.5;
zRes = 1;
vToMPa = 0.8; %500kHz = 0.8V/MPa       2.25MHz = 0.92V/MPa


%reading and processing scan data
M = readmatrix(file);

x = M(1:end,1);
y = M(1:end,2);
z = M(1:end,3);
aV = M(1:end,4);
a1 = aV/vToMPa; %(MPa) convert raw data to MPa
pha = (M(1:end,5) - (max(M(1:end,5))-180)) * (pi/180);

x1 = min(x)-88:xRes:max(x)-88;
y1 = min(y)-160.5:yRes:max(y)-160.5;

[X,Y] = meshgrid(x1,y1);



% rearrange scan data from vector into 2D matrix
ptsPerLayer = length(x1)*length(y1);
numLayers = length(z)/ptsPerLayer;
for layer = 1:numLayers
    for i = 1:length(y1)
        A1(i,:) = a1((i-1)*length(x1)+1 + (ptsPerLayer*(layer-1)):i*length(x1) + (ptsPerLayer*(layer-1)));
        Pha(i,:) = pha((i-1)*length(x1)+1 + (ptsPerLayer*(layer-1)):i*length(x1) + (ptsPerLayer*(layer-1)));
        if mod(i,2) == 0
            A1(i,:) = flip(A1(i,:));
            Pha(i,:) = flip(Pha(i,:));
        end  
    end
    
% plot resulting data    
%     figure
%     surf(X,Y,Pha, 'edgecolor','none')
%     xlabel('X (mm)')
%     ylabel('Y (mm)')
%     axis equal
%     view(2)
    
    figure
%     Pha1 = imgaussfilt(Pha,1);
    surf(X,Y,Pha, 'edgecolor','none')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap jet
    caxis([-pi pi])
    colorbar('Ticks',[-pi,-pi/2,0,pi/2,pi],'TickLabels',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    axis equal
    view(2)

    figure
%     A1 = imgaussfilt(A1,1);
    surf(X,Y,A1, 'edgecolor','none')
    xlabel('X (mm)')
    ylabel('Y (mm)')
    colormap hot
%     caxis([0.006 0.022])
    colorbar
    axis equal
    view(2)
    
%     figure
%     A2 = imgaussfilt(A1,2);
%     surf(X,Y,A2, 'edgecolor','none')
%     xlabel('X (mm)')
%     ylabel('Y (mm)')
%     axis equal
%     view(2)
end





