%% Template for Homework 1 Question 2
clear all
close all
clc

%% Instructions
%  You will first load satpos.mat which is offered on the canvas site. This
%  is a struct that contains 6 fields:
%
%  x: 288x31 matrix that contains the x ECEF coordinates for 31 satellites
%  over 288 epochs
%  y: 288x31 matrix that contains the y ECEF coordinates for 31 satellites
%     over 288 epochs
%  z: 288x31 matrix that contains the z ECEF coordinates for 31 satellites
%     over 288 epochs
%  numSats: The number of active satellites in the constellation
%  numEpochs: The number of epochs (time stamps) in the recorded data. This
%             data was taken at 300s intervals for 24hrs
%  time: This variable is a 288x1 vector that has the time stamp for each
%        epoch. This is given in hours.
%
%  In addition to satpos, your position in lat, long, height and ECEF is 
%  given
%
%  Throughout the code, there are regions that indicate where you should
%  write your own code in.
%
%  Deliverable: There should be one plot for the elevation, 4 plots for the
%  number of satellites in view, and 4 plots for the DOPs.

%% Given information
load('satpos')

mypos_LLH = [37.4241, 237.8339, 0]; % deg, deg, meters
mypos_ECEF = [-2699.958e3, -4293.091e3, 3854.878e3];    % meters

% Extract lat and long and convert to radians
lat = mypos_LLH(1)*pi/180;
long = mypos_LLH(2)*pi/180;

% Compute rotation matrix from ECEF to ENU
R_L = [-sin(long), cos(long), 0;...
    -sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat);...
    cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];

% Compute satellite position in ENU frame
for i = 1:satpos.numSats
    for j = 1:satpos.numEpochs
        [sat_ENU(i).ENU(j,:)] = R_L*[satpos.x(j,i) - mypos_ECEF(1);...
            satpos.y(j,i) - mypos_ECEF(2); satpos.z(j,i) - mypos_ECEF(3)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin Assignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part a) Compute the elevation angle of all of the satellites

% From here, the satellite positions have been converted into the ENU
% frame. Sat_ENU is a struct with 31 fields, one for each satellite. In
% each field there is a 288x3 matrix where the rows correspond to the epoch
% and the columns correspond to E, N, and U, respectively. Familiarize 
% yourself with the structure of sat_ENU in order to continue with the 
% assignment.

% Appendix 4.A in the text is useful for this part of the problem
% It can be found at the end of chapter 4

% Compute elevation angle of satellites vs. time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:satpos.numSats
    for j = 1:satpos.numEpochs
        elevation(i, j) = asin(sat_ENU(i).ENU(j,3)/sqrt(...
         sat_ENU(i).ENU(j,1)^2 + sat_ENU(i).ENU(j,2)^2 + ...
         sat_ENU(i).ENU(j,3)^2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot elevation angles for satellites 1, 10, 20, and 30
% Plot the four elevations on one plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','white')
sats = [1 10 20 30];
for i = 1:length(sats)
    subplot(length(sats), 1, i)
    box on
    plot(elevation(sats(i)))
    grid on
    ylim([-2,2])
    xlim([0 satpos.numEpochs])
    ylabel(['el(' num2str(sats(i)) ') [rad]'])
end
xlabel('time [epochs]')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part b) Compute number of satellites in view using different masking angles

% Masking angles of 10, 15, 20, and 25 degrees
angle = 25; 
mask_angle = angle*pi/180; % Degrees converted to radians

% Count number of satellites in view
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSats_in_view = zeros(1, satpos.numEpochs);
for j = 1:satpos.numEpochs
    for i = 1:satpos.numSats
        if elevation(j) > mask_angle
            numSats_in_view(j) = numSats_in_view(j) + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot number of satelites visible vs. time
% There should be four plots, one for each mask angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
figure('Color','white')
box on
plot(numSats_in_view);
xlim([0 satpos.numEpochs])
ylim([0 inf])
title(['Number of visible satellites above ' num2str(angle) ' deg'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute azimuth vs. time (Given to you)
for i = 1:satpos.numSats
    for j = 1:satpos.numEpochs
        azimuth(j,i) = atan2(sat_ENU(i).ENU(j,1), sat_ENU(i).ENU(j,2));
    end
end

%% Part c) Compute the geometry matrix and calculate the various DOPs using
%  each masking angle from part b)

% Section 6.1.2 in the text is useful for this part of the problem

% Compute geometry (G) matrix at each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = ones(satpos.numSats, 4, satpos.numEpochs);
for j = 1:satpos.numEpochs
   for i = 1:satpos.numSats
      if elevation(i,j) > mask_angle
        G(i, 1:3, j) = sat_ENU(i).ENU(j,:)/norm(sat_ENU(i).ENU(j,:));
      else
        G(i, 1:3, j) = [0 0 0];
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the H matrix at each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDOP = zeros(1, satpos.numEpochs);
HDOP = zeros(1, satpos.numEpochs);
VDOP = zeros(1, satpos.numEpochs);

H = zeros(4, 4, satpos.numEpochs);
for j = 1:satpos.numEpochs 
    H(:, :, j) = inv(G(:, :, j)'*G(:, :, j));
    PDOP(j) = sqrt(H(1, 1, j) + H(2, 2, j) + H(3, 3, j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute PDOP, HDOP, and VDOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gtilde = ones(satpos.numSats, 4, satpos.numEpochs);
for j = 1:satpos.numEpochs
   for i = 1:satpos.numSats
       if elevation(i,j) > mask_angle
           Gtilde(i, 1:3, j) = [cos(elevation(i,j))*sin(azimuth(j,i));...
                           cos(elevation(i,j))*cos(azimuth(j,i));...
                           sin(elevation(i,j))]';
       else
           Gtilde(i, 1:3, j) = [0 0 0];
       end
   end
end

Htilde = zeros(4, 4, satpos.numEpochs);
for j = 1:satpos.numEpochs 
    Htilde(:, :, j) = inv(Gtilde(:, :, j)'*Gtilde(:, :, j));
    HDOP(j) = sqrt(Htilde(1, 1, j) + Htilde(2, 2, j));
    VDOP(j) = sqrt(Htilde(3, 3, j));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot results
% PDOP, HDOP, and VDOP should be on the same plot. There should be 4 total
% DOPs plots for each of the mask angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Your code here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','white')
title(['DOPs for ' num2str(angle) ' deg mask angle'])
box on; grid on; hold on
plot(PDOP)
plot(HDOP,':')
plot(VDOP,'-.')
xlim([0 satpos.numEpochs])
ylim([0 inf])
legend('PDOP','HDOP','VDOP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












