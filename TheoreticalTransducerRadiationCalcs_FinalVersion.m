clear
close all
clc

res = 500;

f = 500e3; %frequency (Hz = 1/s)
rad = 19e-3; %radius of transducer (m)

rho0 = 1000; %equilibrium density = density of water (kg/m^3)
c = 1480; %speed of sound in water (m/s)
L = c/f; %wavelength (m)
w = f*2*pi; %frequency of signal (rad/s)
K = 2*pi/L; %wave number (1/m)
t = 0.1; %linspace(0,0.1,500); %time (s)

volt = 2.8; %Vpp (V)
dh = (250e-12)*volt; %max displacement (m)      400e-12 for 2.25MHz elements      250e-12 for 500kHz elements
U0 = w*dh; %max speed of transducer surface (m/s)

r = linspace(0,rad,65);
dr = r(2)-r(1);
phi = linspace(0,2*pi,65);
dPhi = phi(2)-phi(1);
x = linspace(-00.5,00.5,res);%-0.025,0.025,50);
z = linspace(0,1,res);%0.007,0.057,50);
p = zeros([length(x),length(z),length(t)]); %complex pressure


%Numerical Approximation of Exact Equation
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

[X,Z] = meshgrid(x,z);
X = X.';
Z = Z.';

PlotTitle = sprintf('Transducer Acoustic Field (MPa) \n f = %.2fkHz, Vpp = %.2f', f/(10^3), volt);

figure
H = surf(X,Z,PMPa(:,:,1));
axis equal
colorbar
set(gca,'ColorScale','log','ydir','reverse')
set(H, 'edgecolor', 'none')
xlabel('X (m)')
ylabel('Z (m)')
title(PlotTitle)
view(2)




%Analytic Calcs for comparison
size = 1.00;
pts = 1000;

%Axial Simplification
rAxial = linspace(0,size,80000);
Paxial = 2.*rho0.*c.*U0.*abs(sin(0.5.*K.*rAxial.*(sqrt(1+(rad./rAxial).^2)-1)));

%Far Field Approximation
for i = 1:length(x)
    for k = 1:length(z)
        [theta,R] = cart2pol(z(k),x(i));
        Pax = 0.5*rho0*c*U0*(rad./R)*K*rad;
        v = K*rad*sin(theta);
        H = abs(2*besselj(1,v)/v);
        Pfar(i,k) = Pax*H;
    end
end


PaxialMPa = Paxial./(10^6); %MPa
PfarMPa = Pfar./(10^6); %MPa

figure
plot(rAxial,PaxialMPa)
grid on
title('Axial Pressure from Exact Eq')
xlabel('Radial Distance (m)')
ylabel('Pressure (MPa)')

% figure
% plot(rAxial,PaxialMPa)
% xlim([0,0.2])
% grid on

figure
H = surf(X,Z,PfarMPa(:,:,1));
title('Far Field Approximation (MPa)')
xlabel('X (m)')
ylabel('Z (m)')
colorbar
axis equal
set(gca,'ColorScale','log','ydir','reverse')
set(H, 'edgecolor', 'none')
view(2)

figure
PercentError = (abs(PMPa-PfarMPa)./abs(PMPa))*100;
AreaAvgError = sum(PercentError.*(max(z)/length(z))^2, 'all') / (max(x)*max(z))
PercentError = abs(filloutliers(PercentError,'linear'));
H = surf(X,Z,PercentError);
title('% Error Between Numeric Approx. and Far Field Approx. (MPa)')
xlabel('X (m)')
ylabel('Z (m)')
colorbar
axis equal
set(gca,'ydir','reverse')
set(H, 'edgecolor', 'none')
view(2)


figure
PercentError(PercentError>80) = 80;
H = surf(X,Z,PercentError);
% title('% Error Between Numeric Approx. and Far Field Approx. (MPa)')
xlabel('X (m)')
ylabel('Z (m)')
colorbar
axis equal
set(gca,'ydir','reverse')
set(H, 'edgecolor', 'none')
view(2)







