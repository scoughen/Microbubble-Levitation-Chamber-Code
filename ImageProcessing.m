%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script takes reads an image, detects any circles (bubbles) present,
%  and measures the radius of the circle in both pixels and um.
%
%  The parameters that need to be set are:
%    picName = name of video file (including file extenstion)
%    pxlSize = size of pixel in m
%
%  S. Coughenour - Nov. 17, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

picName = 'Image5.jpg';
pxlSize = 2*10^-6; %m

I = imread(picName);
% imshow(I)

[centers,radii] = imfindcircles(I,[20 60],'ObjectPolarity','dark', ...
          'Sensitivity',0.95,'Method','twostage');
      
figure
imshow(I)
h = viscircles(centers,radii);

umRadius = radii*2; %um



