%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script calculates the modes of a rectangular tank filled with water.
%
%  The parameters that need to be set are:
%    Lx = x length of tank
%    Ly = y length of tank
%    Lz = height of water
%    numModes = number of shape modes
%
%  S. Coughenour - Nov. 17, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
format compact
clc


Lx = 24.9; %cm  ---  Length of tank
Ly = 46.8; %cm  ---  Width of tank
Lz = 50.05-10.5; %cm  ---  Height of the water
numModes = 1000;


c = 1.48*10^5; %cm/s

% Calculating all of the shape modes - shape mode = (l,m,n)
for l = 0:100
    for m = 0:100
        for n = 0:100
            f(l+1,m+1,n+1) = (c/2) * sqrt( (l/Lx)^2 + (m/Ly)^2 + (n/Lz)^2 );  % tank modes equation
        end
    end
end


fvec = reshape(f,1,[]); % convert 3-d matrix into vector
fvecSort = sort(fvec/1000); % sort fvec and convert to kHz
fvecSortTrim = fvecSort(2:numModes); % extract 1000 lowest frequency modes



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot 1000 lowest frequency modes
figure
plot(fvecSortTrim,'.')
ylabel('Resonant Frequency (kHz)')


% Other plotting options

% figure
% hold on
% ylabel('Resonant Frequency (kHz)')
% for i = 1:length(f)
%     plot(f(:,:,i)/1000,'.')
% end
% 
% 
% figure
% hold on
% ylim([497*10^3,503*10^3])
% 
% for i = 1:length(f)
%     plot(f(:,:,i),'.')
% end
% 
% plot(118:length(f),500*10^3*ones(size(118:length(f))))



