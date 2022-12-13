%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script caculates the resonant frequency of a bubble for a given
%  radius.
%
%  The parameters that need to be set are:
%    R0 = the static radius of the bubble
%    h = the distance from the bubble to the surface of the water
% 
%  S. Coughenour - Dec. 1, 2022
%  M. Calvisi - Dec. 1, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

R0 = 106 *10^-6; %static radius (um)
h = 40 *10^-3; %distance of bubble from water surface (mm)

gama = [1,1.4]; %ratio of specific heat of a gas at constant pressure to that at constant volume
rho = 997; %fluid density (kg/m^3)
p0 = rho*9.81*h + 101325; %hydrostatic liquid pressure
w = (1/R0)*sqrt(3*gama*p0/rho); %angular frequency
f = w/(2*pi) %frequency (Hz)