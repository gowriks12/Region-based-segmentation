clc; clear all; close all;
%% "RIGHT CODEEEEE" final final finaaaalllllllllllllll
img = imread("dif00016.pgm"); %Input image
% figure();
% imshow(img);
%% Region-growing
region=zeros(size(img)); %Initializing a matrix to hold region values
regions(1) = img(1,1); %List to hold pixel intensity value for each region and initializing it to first pixel
current_region = 1; %Initializing the current region
for i=1:size(img,1)
    for j=1:size(img,2) 
        found=0; %reset Classiffication flag each iteration
        if(img(i,j) == regions(current_region)) %If the pixel belongs to the current region
            region(i,j) = current_region; %Set the value of current region to the  region pixel
            found=1; %Set classified flag
        else %If pixel does not belong to the region,
            for k=1:size(regions,2) % Check for intensity value that matches the pixel in the 
                if(regions(k)==img(i,j)) %regions list. If found, set the current region to k
                    current_region=k;
                    found=1; %Set classified flag
                end
            end
            if(found==0) %If the regions list does not have that value
                current_region = current_region+1; %Increment the current region value
                regions(current_region)=img(i,j); %Add the intensity value to the list
            end        
            region(i,j) = current_region; %Set the value of current region to the  region pixel
        end
    end
end  
%% 8-connected
neigh = [1 1 1;1 0 1;1 1 1]; %8 surrounding neighbors
for x=1:size(region,1)-1 %Iterate through the region matrix
    for y=1:size(region,2)-1
        if (region(x,y)==2) %If the pixel belongs to object
            res = ((region(x-1:x+1,y-1:y+1).*neigh)==1);%Observe 8 neighbours of current pixel check how many are equal to 1         
            if(sum(res)>=5) %If more than or equal to 5 belong to region 1 or background,
                region(x,y)=1; %Set the region to 1(background)
            end
        end
    end
end

%% Calculating moments
clear M00; clear M01; clear M10; clear M11; clear M02; clear M20;
M00 = 0;M01 = 0;M10 = 0;M11 = 0;M02=0;M20=0;
fimg = 1; %Function f(x,y)=1 for the object
j=1; 
for x= 1:size(img,1) %Iterate through the region matrix
    for y = 1:size(img,2)
        if(region(x,y)==2) %If the region value is 2
            x_coords(1,j)=x; %Copy the x value
            y_coords(1,j)=y; %copy the y value
            j=j+1;
            M00 = M00 + (x^0)*(y^0)*fimg; %Calculate moments
            M01 = M01 + (x^0)*(y^1)*fimg;
            M10 = M10 + (x^1)*(y^0)*fimg;
            M11 = M11 + (x^1)*(y^1)*fimg; 
            M02 = M02 + (x^0)*(y^2)*fimg;
            M20 = M20 + (x^2)*(y^0)*fimg;
        end
    end
end
%% Computations

%Centroid
y_bar = round(M10/M00); 
x_bar = round(M01/M00);

%Orientation
U11 = M11 - ((M10*M01)/M00);
U20 = M20 - ((M10^2)/M00);
U02 = M02 - ((M01^2)/M00);
x_theta = 2*U11/sqrt(4*(U11^2)+(U20-U02)^2);
y_theta = (U20-U02)/sqrt(4*(U11^2)+(U20-U02)^2);
thetad = 0.5*atan2d(x_theta,y_theta);
thetar = pi*thetad/180;
R = [ cos(thetar)   sin(thetar)
     -sin(thetar)   cos(thetar)];
 
%Maximum length 
max_len = x_coords(length(x_coords))-(x_coords(1));


x = [min(y_coords):max(y_coords)]; %Length of vertical line
y_len = [min(x_coords):max(x_coords)]; %Length of horizontal line
%% Plotting
m=R(1,1)/R(1,2); %Slope = cos(theta)/sin(theta)
c = y_bar-m*x_bar; %Intercept of the line passing through centroid
% Vertical Principal axis
verti = m.*x + c; %Equation of vertical line

%Horizontal Principal axis
thetad = thetad+90; %Horizontal axis is perpendicular to vertical line
thetar = pi*thetad/180;
m1 = sin(thetar)/cos(thetar); %Slope in this case is sin(theta)/cos(theta)
% m1 = 1/(m); % Slope of the like that is perpendicular to vertical line

c1 = x_bar-m1*y_bar; %Intercept of the line passing through centroid
hori = (m1.*y_len +c1); %Equation of horizontal line

%% Plotting
figure();
imshow(img); %Show image
hold on
grid on
plot(x_bar, y_bar, 'b*'); %Plotting centroid
plot(x, verti,'r-'); %Plotting vertical line
plot(hori,y_len); %Plotting horizontal line
hold off
grid off
%% Printing the stats
disp("Title");disp("Value");
disp("m00");disp(M00);
disp("m01");disp(M01);
disp("m10");disp(M10);
disp("m11");disp(M11);
disp("m02");disp(M02);
disp("m20");disp(M20);
disp("u11");disp(U11);
disp("u20");disp(U20);
disp("u02");disp(U02);
disp("x bar");disp(x_bar);
disp("y bar");disp(y_bar);
disp("x theta");disp(x_theta);
disp("y theta");disp(y_theta);
disp("orientation(in deg)");disp(thetad);
disp("Slope of vertical line");disp(m);
disp("Slope of horizontal line");disp(m1);
