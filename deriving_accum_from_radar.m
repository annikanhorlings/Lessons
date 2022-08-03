% This script derives snow accumulation from a shallow ice-penetrating
% radar profile collected in Antarctica.

%Author: Annika Horlings
%University of Washington
%Last updated: 15 June 2022

%Directions: Run each cell by selecting the cell then pressing the 
%"Run Section" button at the top of the MATLAB editor. There are questions
%within some cells and these are highlighted with %%% lines. In addition,
%you are encouraged to change some lines of code; these are also clearly
%marked with %%% lines. 


%% 1 - picking the internal layers within the radar data

% (1) Pull up the "iterm" application on this computer.

% (2) Enter this line of code to change directories and locate the radar
    % data (without the % at the beginning).
    % cd /Users/annikahorlings/Documents/TeachingOutreach_and_Coordinating/JIRP2022/Workshops/Snow_accumulation_from_radar

% (3) Load the data into ImpDAR, a processing and interpretation
    % toolbox for ice-penetrating radar data (Lilien et al., 2019). Enter
    % this line of code into the "iterm" command line (without the % at the 
    % beginning):
    % imppick shallow_radar_data.mat
    
% (4) You should see now a new window with the radar data displayed. Locate 
    % the left-hand panel and change the frequency to 100 and pick 
    % number to 1.

% (5) How do we start picking a layer? Select "New Pick" and switch the
    % mode from "Select Mode" to "Edit Mode". Now you're ready to pick the
    % layer! Locate the internal layer that is exceptionally bright and
    % that is just greater than 1.00 microseconds in two-way travel time.
    % This is the layer you want to pick across the profile.
    
% (6) Once you have picked the layer, close the window and save the data
    % with a new name in the same folder.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is an internal layer represent within the radar data? Why does it
% appear? Why are some layers brighter than others?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 - Calculating the accumulation rate 

% Overview - How do we calculate the accumulation rate? 
% (a) Convert the two-way travel time to depth. This requires knowing the
    % velocity of the radar wave in the firn or ice.
    
% (b) Determine the age of the internal layer that is picked (because, of
    % course, the accumulation rate is a RATE ...). 

% (c) Use the equation below to derive the time-averaged accumulation rate:
    % accumulation = (depth/age)*(mean density)/(density of ice));

%Let's get started!

%% 2.1 - Convert the two-way travel time of the layer to depth
cd('/Users/annikahorlings/Documents/TeachingOutreach_and_Coordinating/JIRP2022/Workshops/Snow_accumulation_from_radar');
p = load('shallow_radar_data_example_pick.mat'); picksin = p.picks.samp2;
datain = p.data; tn = p.trace_num; TWTin = p.travel_time;
[idxcrop, TWT, data] = adjustTWT(TWTin, datain);
[idx, picks] = adjustpicks(picksin, TWTin, datain);
distance = zeros(1, length(p.x_coord));
for i = 2:length(p.x_coord)
    distance(1) = 0;
    distance(i) = (sqrt((p.x_coord(i) - p.x_coord(i-1)).^2 + ...
        (p.y_coord(i) - p.y_coord(i-1)).^2))./1000 + distance(i-1);
end

%find the TWTT that corresponds with the picks
[m n] = size(picks);
for k = 1:m
    for l = 1:n
        if isnan(picks(k, l)) == 1
            TWTT(k,l) = NaN;
        else
            TWTT(k,l) = TWT(picks(k, l));
        end
    end
end
TWTT_sec = TWTT./(10^6); %convert to seconds

for i = 1:length(TWTT)
       depth(i) = TZ_firn_bdot_T(TWTT_sec(i), 0.12, -40);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The last few lines of code use a function that converts the two-way
% travel time to depth using a firn depth-density profile defined by the
% last two entries (anticipated snow acccumulation rate and surface
% temperature, respectively) to get the radar velocities.
% You may realize that this seems circular and we'll be looking into this 
% later.

% Question - How deep is the internal layer that you traced?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure(1);
plot(distance, depth, 'LineWidth', 2);
xlabel('distance (km)');
ylabel('depth (m)');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Reverse');
f1.Position = [217   524   912   408];


%% 2.2 Determine the age of the internal layer that you traced

% This is the depth-age scale that we have:
agechem_hd = table2array(readtable('US_ITASE-02-4_2013.csv'));
depthcore = agechem_hd(:, 1);
agecore = agechem_hd(:, 2);
depthcoreextended = vertcat(depthcore, linspace(72.05, 75, 500)');
agecoreextended = interp1(depthcore, agecore, depthcoreextended, ...
    'linear', 'extrap');
meandepth = nanmean(depth);
[minvalue, idxdepth] = min(abs(meandepth - depthcoreextended));
age = agecoreextended(idxdepth);
ageold = 2019 - age;

figure(2);
plot(agecoreextended, depthcoreextended, 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Reverse');
xlabel('Date');
ylabel('Depth (m)');

%% 2.3 We now need to know how much mass is above the layer
[z, rho] = HL_analytic_adj((273+(-40)), 0.12, 0.35);
zdiscr = (z(2, 1) - z(1, 1)); %Calculate the ice above the depth of internal layers
for j = 1:m
     [valz(j, 1), idxz(j, 1)] = min(abs(z - meandepth(j)));
for i = 1:idxz(j) %discretize and sum up the amount of ice above this depth
    m_parcel(j, i) = rho(i, 1)*zdiscr*1*1; %calculate the mass in each discretized parcel
end
    m_all(j) = sum(m_parcel(j, :)); %calculate all ice mass above given depth
    rho_mean(j) = m_all(j)/z(idxz(j))/1/1; %divide by the total ice thickness to get a mean density of the column
end

%% 2.4 Now compute the accumulation rates
accum = depth./ageold.*(rho_mean/0.917);

f3 = figure(3);
plot(distance, accum, 'LineWidth', 2, 'Color', [204/255 102/255 0]);
xlabel('distance (km)');
ylabel('accumulation rate (m yr^{-1})');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Reverse');
f3.Position = [217   524   912   408];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Why do you think there's so much variation in th accumulation rate along 
% this profile? What are some of the many factors that influence snow
% accumulation at a specific location?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3 - Thinking about processes

% Let's think about two processes: local topographic impacts on wind
% redistribution and the broader-scale influence of orography. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What is the orographic effect? Also, what is one way we could think about
% local topographic influence, i.e., how can we quantitatively compute the
% convexity or concavity of the topography? (Think calculus here.) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3.1 Local toporaphic influence 

f4 = figure(4);
subplot(3, 1, 1);
plot(distance, smooth_ndh(p.elev, 50, 0), 'LineWidth', 2, 'Color', [204/255 102/255 0]);
xlabel('distance (km)');
ylabel('surface elevation (m)');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Normal');

subplot(3, 1, 2);
plot(distance, smooth_ndh(accum, 5, 0), 'LineWidth', 2, 'Color', [102/255 102/255 0]);
xlabel('distance (km)');
ylabel('accumulation rate (m yr^{-1})');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Normal');

subplot(3, 1, 3);
plot(distance, smooth_ndh(gradient(gradient(smooth_ndh(p.elev, 50, 0))),...
    200, 0), 'LineWidth', 2, 'Color', [102/255 50/255 0]);
xlabel('distance (km)');
ylabel('concavity/convexity');
grid on;
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Ydir', 'Normal');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What do you think? Is there a clear correlation? What else must we
% consider when doing this analysis?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Notes: "Accumulation rates are higher in topographic depressions. Over the 
% Antarctic ice sheet, King and others (2004) showed that topography, 
% snow precipi- tation, wind direction and speed resulted in significant 
% snow redistribution that leads to high spatial variability in snow 
% accumulation over distances as short as a few kilometers. The degree to 
% which surface slope and accumulation rate correlate depends on how the 
% GPR profile is oriented relative to the direction of the dominant wind 
% ?The second derivative of the surface-elevation profile
% yields the topographic curvature corresponding to the convexity or 
% concavity of the terrain. Positive values represent concavity, where the 
% terrain has a local depres- sion, and negative values represent convexity, 
% where the terrain has a local high." (Miege et al., 2013)

%% 3.2 Orographic influence
% Let's look at the broader-scale toporaphy of the region (%this part of
% the code can take a little time to run).

cd('/Users/annikahorlings/Documents/Projects/Layer_thinning/Paper/Final_RunsMARCH2020/Horizontal_strain_rates_maps/rema');
[xr,yr,Zr] = rema_data('xy', [-480*1000 -340*1000], [-170*1000 -80*1000]);
xr; % xim: the x axis of the underlying image
yr; % yim: the y axis of the underlying image
Zr; %get(gcf, 'Colormap'); % im: the n x m x 1 values for the ArcticDEM hillshade
minc = 2520;% minc/maxc: the limits for black and white in your hillshade
maxc = 2650;
im_rgb = colorlock_mat(Zr,gray,[minc maxc]);
corelat = -86.5025;
corelong = -107.990313;
[corex, corey] = polarstereo_fwd(corelat, corelong, 6378273, 0.081816153, -70, 0);

f5 = figure(5);
imagesc(xr/1000, yr/1000, Zr);
hold on;
[C, h] = contour(xr/1000, yr/1000, Zr, ...
    [2480 2500 2520 2540 2560 2580 2600 2620 2640], 'Color', ...
    'k', 'LineWidth', 0.25, 'LineStyle', '-');
set(gca, 'Ydir', 'Normal');
hold on;
xlim([-480 -340]);
ylim([-170 -80]);
hold on;
hold on;
xlabel('PS E (km)');
ylabel('PS N (km)');
colormap(gray);
colorbar;
colr.Label.String = 'Accumulation rate (m i.e. per yr)';
caxis([2480 2650]);
set(gca, 'LineWidth', 3)
set(gca, 'FontSize', 12);
hold on;
plot(p.x_coord/1000, p.y_coord/1000, 'Color', [0 102/255 204/255],...
    'LineWidth', 2);
hold on;
plot(corex/1000, corey/1000, 'p', 'LineWidth', 2, 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'w', 'MarkerSize', 14);

axes('Position',[.13 .13 .2 .2]);
imagesc(xr, yr, Zr);
set(gca, 'Ydir', 'Normal');
caxis([-100 -50])
antbounds('gl','black', 'LineWidth', 1.5)
antbounds('coast','black', 'LineWidth', 2);
hold on;
%antbounds('shelves','Color', [0.5 0.5 0.5]);
hold on;
plot(corex/1000, corey/1000, 'd', 'LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [153/255 76/255 0], 'MarkerSize', 8);
box off;
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you knew the winds for most storm came from the north-northwest, how 
% do you think orography would affect the snow accumulation of this region
% and/or of this profile? (Blue line is the radar profile location.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








