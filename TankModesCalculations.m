clear
close all
format compact
clc

Lx = 24.9; %cm
Ly = 46.8; %cm
Lz = 50.05-10.5; %cm ---------- Should be water height!!!!!!!!!!!!
c = 1.48*10^5; %cm/s

for l = 1:118
    for m = 1:118
        for n = 1:118
            f(l,m,n) = (c/2) * sqrt( (l/Lx)^2 + (m/Ly)^2 + (n/Lz)^2 );
        end
    end
end

fvec = reshape(f,1,[]);
% fvecIndex = 
fvecSort = sort(fvec/1000); %outputs in kHz


figure
plot(fvecSort,'.')
ylabel('Resonant Frequency (kHz)')

figure
hold on
ylabel('Resonant Frequency (kHz)')
for i = 1:length(f)
    plot(f(:,:,i)/1000,'.')
end


figure
hold on
ylim([497*10^3,503*10^3])

for i = 1:length(f)
    plot(f(:,:,i),'.')
end

plot(118:length(f),500*10^3*ones(size(118:length(f))))



