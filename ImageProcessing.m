clear
close all
clc

I = imread('Image5.jpg');
% imshow(I)

[centers,radii] = imfindcircles(I,[20 60],'ObjectPolarity','dark', ...
          'Sensitivity',0.95,'Method','twostage');
      
figure
imshow(I)
h = viscircles(centers,radii);

umRadius = radii*2; %um