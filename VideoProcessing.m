clear
close all
clc

vid = VideoReader('BubbleOscillation3_27.02kHz_100Vpp_80000fps.mp4'); %might also work 'video.avi'
numFrames = vid.NumberOfFrames;
fps = 80000;
dt = 1/fps;
t = 0;

freq = 27.02*10^3; %Hz
period = 1/freq;
period10 = 10*period;
endFrame = ceil(period10/dt);

for i = 1:1:endFrame%100
    frames = read(vid,i);
%     imwrite(frames,['Image' int2str(i), '.jpg']);
%     im=image(frames);
%     name = sprintf("Image%d.jpg", i);
    I = frames; %imread(name);
    [center,radius] = imfindcircles(I,[20 41],'ObjectPolarity','dark','Sensitivity',0.88,'Method','PhaseCode');
    radii(i) = radius;
    centers(i,1:2) = center;
    
    tt(i) = t;
    t = t+dt;
    
%     figure
%     hold on
%     imshow(I)
%     h = viscircles(center,radius);
end

radii_um = radii*2*10^-6;

figure
plot(tt,radii_um)
xlabel('Time (s)')
ylabel('Radius (m)')


%Curve Fitting
% FT = fittype('abs(1.9400 * sin(f*(x+C))) + 24.8203'); %fittype('abs(A*sin(f*(x+C)))+B');  % y = f(x)
% f0 = fit(tt',radii',FT)
% 
% hold on
% plot(f0)
% 
% 
% customFit = abs( (max(radii)-min(radii)) * sin(freq*2*pi * (linspace(0,max(tt),1000) - 0.45e-5))) + min(radii);
% plot(linspace(0,max(tt),1000),customFit)





