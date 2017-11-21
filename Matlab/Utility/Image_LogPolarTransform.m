function retval_imgLogPolar2D = Image_LogPolarTransform(imgCartesian2D,num_angles)

imgCartesian2D(isnan(imgCartesian2D))=0; 


x0=size(imgCartesian2D,1)/2;
y0=size(imgCartesian2D,2)/2;
[y_coordinates, x_coordinates, intensity]=find(imgCartesian2D);
x_coordinates=x_coordinates-x0;
y_coordinates=y_coordinates-y0;
theta  = atan2(y_coordinates,x_coordinates);
radius = sqrt(x_coordinates.^2+y_coordinates.^2);

% Determine the minimum and the maximum x and y values:

rthresh=min(x0-0,y0-0);

indices_inside=find(radius<rthresh);
theta  = theta (indices_inside);
radius = radius(indices_inside);
x_coordinates=x_coordinates(indices_inside);
y_coordinates=y_coordinates(indices_inside);
intensity = intensity(indices_inside);

tmin = min(theta); tmax = max(theta);
rmin = min(radius); rmax = max(radius);


% Define the resolution of the grid:
rres=rmax; % # of grid points for R coordinate
tres=num_angles; % # of grid points for theta coordinate

F = TriScatteredInterp(radius,theta,intensity,'natural');

%Evaluate the interpolant at the locations (rhoi, thetai).
%The corresponding value at these locations is Zinterp:

[radius_i,theta_i] = meshgrid(linspace(rmin,rmax,rres),linspace(tmin,tmax,tres));
intensity_i = F(radius_i,theta_i);


retval_imgLogPolar2D = intensity_i;


end

% function imC = Polar2Im(imP,W,method) 
% %Polar2Im turns a polar image (imP) into a cartesian image (imC) of width W 
% %method can be: '*linear', '*cubic', '*spline', or '*nearest'. 
% imP(isnan(imP))=0; 
% w = round(W/2); 
% xy = (1:W-w); 
% [M N P]= size(imP); 
% [x y] = meshgrid(xy,xy); 
% n = round(N/4); 
% rr = linspace(1,w,M); 
% W1 = w:-1:1; 
% PM = [2 1 3;1 2 3;2 1 3;1 2 3]; 
% W2 = w+1:2*w; 
% nn = [1:n; n+1:2*n; 2*n+1:3*n; 3*n+1:N;]; 
% w1 = [W1;W2;W2;W1]; 
% w2 = [W2;W2;W1;W1]; 
% aa = linspace(0,90*pi/180,n); 
% r = sqrt(x.^2 + y.^2); 
% a = atan2(y,x); 
% imC= zeros(W,W,P); 
% for i=1:4 %turn each quarter into a cartesian image 
% imC(w1(i,:),w2(i,:),:)=permute(interp2(rr,aa,imP(:,nn(i,:))',r,a,method),PM(i,:)); 
% end 
% imC(isnan(imC))=0;