%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script takes reads an image, extracts the r(theta) radial
%  description of the bubble, and calculates the nonspherical shape modes.
%
%  The parameters that need to be set are:
%    picName = name of image file (including file extenstion)
% 
%  S. Coughenour - Dec. 1, 2022
%  M. Calvisi - Dec. 1, 2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%Load image
%%%%%%%%%%%
picName = 'NonsphereImage.jpg';
I = imread(picName);
imshow(I)


%Prepare image for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = rgb2gray(I); %convert image to gray scale
BWim = imcomplement(imbinarize(I)); %convert image to binary and invert
% figure  %optional display of current picture
% imshow(BWim)
BWim1 = bwareafilt(BWim,1); %remove all but the largest region
figure  %optional display of current picture
imshow(BWim1)


%Find the initial centroid of the bubble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = regionprops(BWim1,'centroid'); %find centroid
centroids1 = cat(1,s.Centroid); %convert format from struct to double vector
centroid1 = [centroids1(1,1),centroids1(1,2)]; %extract only first centroid in case multiple were found
figure  %show the processed image of the bubble and mark its centroid
imshow(BWim1)
hold on
plot(centroids1(1,1),centroids1(1,2),'b*')


%Find approx axis of symmetry
%%%%%%%%%%%%%%%%%%%%%%%
o = regionprops(BWim1,'orientation'); %find symmetry line
orientation = o.Orientation; %convert format from struct to double vector
line([centroid1(1) centroid1(1)+60],[centroid1(2) tand(orientation)*(centroid1(1)+60)]) %show approx axis of symmetry


%rotate image to put axis of symmetry on x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tform = affine2d([cosd(orientation) -sind(orientation) 0; sind(orientation) cosd(orientation) 0; 0 0 1]); %create rotation matrix
BWimRot = imwarp(BWim1,tform); %rotate image using rotation matrix
figure %show now rotated image
imshow(BWimRot)


%Find the centroid of the bubble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = regionprops(BWimRot,'centroid'); %find centroid
centroids = cat(1,s.Centroid); %convert format from struct to double vector
centroid = [centroids(1,1),centroids(1,2)]; %extract only first centroid in case multiple were found
figure  %show the processed image of the bubble and mark its centroid
imshow(BWimRot)
hold on
plot(centroids(1,1),centroids(1,2),'b*')



%finding r(theta)
%%%%%%%%%%%%%%%%%
dTheta = 0.01; %set angle step size
dr = 0.5; %set radius step size
theta = 0:dTheta:pi; %create angle(theta) vector
x = round(centroid(1)); %set starting x and y coord's for search
y = round(centroid(2));
flag = 1; %initialize flag for finding edge of bubble
r = zeros(1,length(theta)); %prealocate r vector

for i = 1:length(theta) %loop through each angle
    while flag == 1 %loop until edge of bubble detected
        r(i) = sqrt(round(x-centroid(1))^2 + round(y-centroid(2))^2); %update radius
        
        dx = dr*cos(theta(i)); %calculate x increment from raduis step size and current angle
        dy = dr*sin(theta(i)); %calculate y increment from raduis step size and current angle

        x = x+dx; %increment x
        y = y+dy; %increment y

        flag = BWimRot(round(y),round(x)); %update flag
    end
    x = round(centroid(1)); %reset x and y to starting coord's
    y = round(centroid(2));
    flag = 1; %reset flag
end

figure %plop extracted edge of bubble
polar(theta,r)
hold on
polar(-theta,r)


% smooth r
%%%%%%%%%%
rSmooth = smoothdata(r); %remove most noise induced by pixilation

figure %plot now smoothed bubble edge
polar(theta,rSmooth)
hold on
polar(-theta,rSmooth)


% Extract R_0
%%%%%%%%%%%%%
Rmin = min(rSmooth); %find minimum radius (inscribing cirle)
Rmax = max(rSmooth); %find maximum radius (circumscribing cirle)
R_0 = mean([Rmin Rmax]); %use average of min and max radii as approximation of sphrical shape mode

figure(5) %plot inscribing and circumscribing circles on bubble image
hold on
viscircles(centroid,Rmin)
viscircles(centroid,Rmax)


% Extract (compute) the shape mode amplitudes from r_s by doing numerical integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 5; %set highest shape mode to calculate
a_n_comp = zeros(N-1,1); %initialze computed shape mode amplitudes to zero
r_s2 = rSmooth - R_0; %subtract R_0 from r_s before computing shape mode amplitudes

for k = 2 : N
    for j = 1 : length(theta)
        a_n_comp(k-1) = a_n_comp(k-1) + r_s2(j) * legendreP( k, cos(theta(j)) ) * sin(theta(j)) * dTheta ;
    end
    a_n_comp(k-1) = ((2*k + 1)/2) * a_n_comp(k-1);
end


r_s = R_0 * ones( 1, length(theta) ); %initialize r_s to be a 1 x length(theta) vector with each component equal to R_0
for i = 2 : N
    r_s = r_s + a_n_comp(i-1) * legendreP(i,cos(theta));
end 

figure %plot bubble shape given by extracted shape modes
polarplot( theta, r_s )
hold on
polarplot( -theta, r_s )

%rDiff = rSmooth-r_s; %calculate difference between measured radius (smoothed) and radius resulting from extracted shape modes
%figure %plot difference between radii
%plot(theta,rDiff)


%Below code added by MLC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rDiff_norm = ( rSmooth - r_s ) / R_0; %calculate NORMALIZED difference between measured radius (smoothed) and radius resulting from extracted shape modes
figure %plot NORMALIZED difference between radii
plot(theta,rDiff_norm)
grid on
xlabel('Theta (radians)')
ylabel('Normalized Difference')
title('Normalized Difference Between Measured & Computed Radii')

% Plot the NORMALIZED computed (extracted) values of shape mode amplitudes.
a_n_comp_norm = a_n_comp / R_0;     % Normalize the shape mode coefficients with the radius, R_0
modenum = [ 2 : N ];                % Vector of mode numbers (P_1 mode is neglected)
figure
plot( modenum, a_n_comp_norm, 'bo' )
grid on
axis( [0 N+1 min(a_n_comp_norm)-0.02 max(a_n_comp_norm)+0.02 ] )      % Scale axes of plot
xlabel('Mode Number (n)')
ylabel('Shape Mode Amplitude')
title('Normalized Shape Mode Ampltudes')

