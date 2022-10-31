clear
close all
clc


%%%%%%%%%%%%%%%%%%%%%%Scan Stuff%%%%%%%%%%%%%%%%%%%%%

file = "SimpleTransducer500kHzHighRes10SampleAvg0SecDelay7mmOffsetFromTransducerNewAmpTry2.csv"; %SimpleTransducer500kHzLowAmplitude10SampleAvg0SecDelay7mmOffsetFromTransducerWithDamperInTank.csv";       
xRes = 0.5;
yRes = 1;
zRes = 0.5;
vToMPa = 0.8; %500kHz = 0.8V/MPa       2.25MHz = 0.92V/MPa

M = readmatrix(file);

x = M([1:end],1);
y = M([1:end],2);
z = M([1:end],3);
aV = M([1:end],4); %(V)
a = aV/vToMPa; %(MPa)
a = [a;0.01;0.01;0.01];

x1 = min(x):xRes:max(x);
% x1 = x1(1:end-1);
y1 = min(y):yRes:max(y);
z1 = min(z):zRes:max(z);

centerX = 83.5;
centerY = 162;
topZ = max(z);

%  Y  Fig 11 Fig 13	Fig 15	AreaAvgError	AreaAvgError1	AreaAvgErrorNorm
% 159	570	  540	280     	271.6           271.9       	86.5
% 160	470	  480	180         277.8           277.4           86.7
% 161	450	  530	200     	287.5           284.4           90.3
% 162	550	  440	180     	290.3           289.2           89.8
% 163	460	  460	210     	294.2           292.1           93.3
% 164	470	  450	220     	295.7           292.4           94.2
% 165	450	  420	220     	289.8           290.5           94.5
% 166	450	  390	200     	290.3           286.8           94.7
% 167	490	  390	380     	284.1           281.8           91.4
% 168	530	  400	200     	274.3           276.1           86.6


x1 = (x1-centerX)/1000;
y1 = (y1-centerY)/1000;
z1 = -(z1-topZ-7)/1000;

xMax = find(x1==-min(x1));
xCrop = x1(1:xMax);

[X,Y,Z] = meshgrid(xCrop,y1,z1);

ptsPerLayer = length(x1)*length(y1);
numLayers = length(z)/ptsPerLayer;

startIndex = 1;

for layer = 1:numLayers
    for i = 1:length(y1)
        endIndex = startIndex+length(x1)-1;
        A(i,:,layer) = a(startIndex:endIndex);
        startIndex = endIndex+1;
        if mod(i,2) == 0
            A(i,:,layer) = flip(A(i,:,layer));
        end  
    end
    A1 = imgaussfilt(A,1);

%     A2 = imgaussfilt(A,2);

end

Anorm = A./max(max(A));

A = flip(A,3);
A1 = flip(A1,3);
% A2 = flip(A2,3);
Anorm = flip(Anorm,3);

ACrop = A(:,1:length(xCrop),:);
A1Crop = A1(:,1:length(xCrop),:);
% A2Crop = A2(:,1:length(xCrop),:);
AnormCrop = Anorm(:,1:length(xCrop),:);

yMidIndex = find(y1==0);

xx = reshape(X(yMidIndex,:,:),length(xCrop),[]);
zz = reshape(Z(yMidIndex,:,:),length(xCrop),[]);
AA = reshape(ACrop(yMidIndex,:,:),length(xCrop),[]);
AA1 = reshape(A1Crop(yMidIndex,:,:),length(xCrop),[]);
% AA2 = reshape(A2Crop(yMidIndex,:,:),length(xCrop),[]);
AAnorm = reshape(AnormCrop(yMidIndex,:,:),length(xCrop),[]);


% Scan Plot
figure
H = surf(xx,zz,AA);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')
set(gca,'ColorScale','linear','ydir','reverse')
set(H, 'edgecolor', 'none')
view(0,90)

% Scan Plot Axis
xMidIndex = find(x1==0);
figure
plot(z1,AA(xMidIndex,:))
xlabel('Z (m)')
ylabel('Amplitude (MPa)')

TopAvgScan = sum(AA(:,1), 'all') / length(A(:,1))
AreaAvgScan = sum(AA, 'all') / (length(A(:,1))*length(A(1,:)))


% % Scan Plot Smoothed
% figure
% H = surf(xx,zz,AA1);
% axis equal
% colorbar
% xlabel('X (m)')
% ylabel('Z (m)')
% set(gca,'ColorScale','linear','ydir','reverse')
% set(H, 'edgecolor', 'none')
% view(0,90)
% 
% % Scan Plot Smoothed Axis
% figure
% plot(z1,AA1(xMidIndex,:))
% xlabel('Z (m)')
% ylabel('Amplitude (MPa)')


% Scan Plot More Smoothed
% figure
% H = surf(xx,zz,AA2);
% axis equal
% colorbar
% xlabel('X (m)')
% ylabel('Z (m)')
% set(gca,'ColorScale','linear','ydir','reverse')
% set(H, 'edgecolor', 'none')
% view(0,90)

% Scan Plot More Smoothed Axis
% figure
% plot(z1,AA2(xMidIndex,:))
% xlabel('Z (m)')
% ylabel('Amplitude (MPa)')


% Scan Plot Normalized
figure
H = surf(xx,zz,AAnorm);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')
set(gca,'ColorScale','linear','ydir','reverse')
set(H, 'edgecolor', 'none')
view(0,90)

% Scan Plot Normalized Axis
figure
plot(z1,AAnorm(xMidIndex,:))
xlabel('Z (m)')
ylabel('Amplitude (unitless)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%Theory Calcs%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 500e3; %frequency (Hz = 1/s)
rad = 19e-3; %radius of transducer (m)

rho0 = 1000; %equilibrium density = density of water (kg/m^3)
c = 1480; %speed of sound in water (m/s)
L = c/f; %wavelength (m)
w = f*2*pi; %frequency of signal (rad/s)
K = 2*pi/L; %wave number (1/m)
t = 0.1; %linspace(0,0.1,500); %time (s)

volt = 9.6; %0.665; %Vpp (V)
dh = (250e-12)*volt; %max displacement (m)      400e-12 for 2.25MHz elements      250e-12 for 500kHz elements
U0 = w*dh; %max speed of transducer surface (m/s)

xMin = min(xCrop);
xMax = max(xCrop);
zMin = min(z1);
zMax = max(z1);

r = linspace(0,rad,65);
dr = r(2)-r(1);
phi = linspace(0,2*pi,65);
dPhi = phi(2)-phi(1);
x = linspace(xMin,xMax,length(xCrop));
z = linspace(zMin,zMax,length(z1));
p = zeros([length(x),length(z),length(t)]); %complex pressure

for i = 1:length(x) %spacial grid
    for k = 1:length(z)
        
        int = 0;
        for m = 1:length(r)                 %integrate 1/rprime*exp(j*(w*t-k*rprime)) over S
            for n = 1:length(phi)
                dA = r(m)*dr*dPhi;
                r1 = r(m);
                r2 = x(i);
                rprime = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(phi(n)-0) + z(k)^2);
                int = int + (1/rprime)*exp(1i*(w*t-K*rprime)) * dA;
            end
        end
        
        p(i,k,:) = 1i*rho0*c*U0/L*int;
        
    end
end

P = abs(p); %Pa
PMPa = P./(10^6); %MPa
Pnorm = PMPa./max(max(PMPa));

[X,Z] = meshgrid(x,z);
X = X.';
Z = Z.';

% Theory Plot
figure
H = surf(X,Z,PMPa(:,:,1));
axis equal
colorbar
set(gca,'ColorScale','linear','ydir','reverse')
set(H, 'edgecolor', 'none')
xlabel('X (m)')
ylabel('Z (m)')
view(2)

% Theory Plot Axis
figure
plot(z,PMPa(xMidIndex,:))
xlabel('Z (m)')
ylabel('Amplitude (MPa)')


% Normalized Theory Plot
figure
H = surf(X,Z,Pnorm(:,:,1));
axis equal
colorbar
set(gca,'ColorScale','linear','ydir','reverse')
set(H, 'edgecolor', 'none')
xlabel('X (m)')
ylabel('Z (m)')
view(2)

% Normalized Theory Plot Axis
figure
plot(z,Pnorm(xMidIndex,:))
xlabel('Z (m)')
ylabel('Amplitude (normalized)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Data Comparison%%%%%%%%%%%%%%%%%%%%
AA = rot90(rot90(AA));
AA1 = rot90(rot90(AA1));
% AA2 = rot90(rot90(AA2));
AAnorm = rot90(rot90(AAnorm));

Diff = abs((AA-PMPa)./PMPa)*100;
Diff = filloutliers(Diff,'linear');
AreaAvgError = sum(Diff.*(max(z)/length(z))^2, 'all') / (max(x)*max(z))

Diff1 = abs((AA1-PMPa)./PMPa)*100;
Diff1 = filloutliers(Diff1,'linear');
AreaAvgError1 = sum(Diff1.*(max(z)/length(z))^2, 'all') / (max(x)*max(z))

% Diff2 = abs((AA2-PMPa)./PMPa)*100;
% Diff2 = filloutliers(Diff2,'linear');
% AreaAvgError2 = sum(Diff2.*(max(z)/length(z))^2, 'all') / (max(x)*max(z))

Diffnorm = abs((AAnorm-Pnorm)./Pnorm)*100;
Diffnorm = filloutliers(Diffnorm,'linear');
AreaAvgErrorNorm = sum(Diffnorm.*(max(z)/length(z))^2, 'all') / (max(x)*max(z))


% Scan Error Plot
figure
H = surf(xx,zz,Diff);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')
set(gca,'ColorScale','linear')
set(H, 'edgecolor', 'none')
view(0,90)

% Scan Error Plot Axis
figure
plot(z,Diff(xMidIndex,:))
xlabel('Z (m)')
ylabel('% Error')


% Scan Error Plot Smoothed
figure
H = surf(xx,zz,Diff1);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')
set(gca,'ColorScale','linear')
set(H, 'edgecolor', 'none')
view(0,90)

% Scan Error Plot Smoothed Axis
figure
plot(z,Diff1(xMidIndex,:))
xlabel('Z (m)')
ylabel('% Error')


% Scan Error Plot More Smoothed
% figure
% H = surf(xx,zz,Diff2);
% axis equal
% colorbar
% xlabel('X (m)')
% ylabel('Z (m)')
% set(gca,'ColorScale','linear')
% set(H, 'edgecolor', 'none')
% view(0,90)
% 
% Scan Error Plot More Smoothed Axis
% figure
% plot(z,Diff2(xMidIndex,:))
% xlabel('Z (m)')
% ylabel('% Error')


% Scan Error Plot Normalized
figure
H = surf(xx,zz,Diffnorm);
axis equal
colorbar
xlabel('X (m)')
ylabel('Z (m)')
set(gca,'ColorScale','linear')
set(H, 'edgecolor', 'none')
view(0,90)

% Scan Error Plot Normalized Axis
figure
plot(z,Diffnorm(xMidIndex,:))
xlabel('Z (m)')
ylabel('% Error')






