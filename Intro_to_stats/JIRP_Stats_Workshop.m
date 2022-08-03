%% Statistics Workshop
%Author: Annika Horlings
%July 2022

% Instructions: Run each sell by clicking on the section and selecting the 
% "Run Section" button in the MATLAB editor. The script is designed so that
% each section will either output a figure illustrating the concept, or
% some text in the command window that will be helpful.

%% 1 Load in data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We have two datasets that we'll be using: (1) the World Glacier Inventory
%and (2) the SUMUP firn density database. Let's check these out and perform
%different exploratory visualizations and statistics on them!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%World Glacier Inventory:
cd('/Users/annikahorlings/Documents/TeachingOutreach_and_Coordinating/JIRP2022/Workshops/Intro_to_stats/');
wgitab = readtable('glacier_inventory_Alaska.csv');
akarea = table2array(wgitab(:, 5)); akelev = table2array(wgitab(:, 20));
akabllen = table2array(wgitab(:, 12)); 

%Firn from SUMUP database:
filename = 'sumup_density_2020.nc'; %file name
ncdisp(filename); %display contents of the nc file
densityin = ncread(filename,'Density'); %read in the density field
depthin = ncread(filename,'Midpoint'); %read in the midpoint depth field
depthin(depthin < 0) = NaN; depth = depthin(~isnan(depthin)); 
densityin(densityin < 0) = NaN; density = densityin(~isnan(depthin)); 
depth = depth(~isnan(density));
density = density(~isnan(density));

%% 2 Data visualizations

%% 2.1 Histograms
% A histogram provides a quick visual insight into how a data set is distributed. 
% The range of possible values is divided into intervals, or bins. Then a bar chart 
% is created, where the height of each bar corresponds to how frequently values 
% in that bin appear in the data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%You can change this number to alter the number of bins:
%Hint: try something between 40 and 200 for this particular dataset.
nbins = 40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
h = histogram(akelev, nbins);
ylabel('Number of Alaskan Glaciers');
xlabel('Mean Elevation of Accumulation Zone (m)');
set(gca, 'FontSize', 14, 'LineWidth', 2);
title('Histogram');
grid on;

%% 2.2 Box Plots
% A box plot is another way to visualize the distribution of a data set. 
% The central box represents the middle 50% of observations, with the red 
% line at the median. The "whisker" lines show the extent of 99% of the 
% data. Remaining outliers are shown individually with red crosses.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% You can try taking only sections of the data here and see what it does 
% to the boxplot:
% Hint: Try taking the 1:1000, 1000:2000, 1000:14288, or 1:14288 (like
% above)
section = 1:1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
b = boxplot(akelev(section));
xlabel('Alaskan Glaciers');
ylabel('Mean Elevation of Accumulation Zone (m)');
set(gca, 'FontSize', 14, 'LineWidth', 2);
title('Box Plot');


%% Scatter Plots
% A scatter plot explores how two variables are related to each other. 
% You can use the scatter function or plot function to create a scatter 
% plot.
figure(3);
s1 = scatter(akelev, akarea, 20, [204/255 102/255 0]);
xlabel('Mean Elevation of Accumulation Zone (m)')
ylabel('Area of Alaskan Glaciers (km^{3})');
ylim([0 500]);
xlim([-10 4000]);
set(gca, 'FontSize', 14, 'LineWidth', 2);
grid on; box on;
title('Scatter Plot #1 - Alaskan Glaciers');
grid on; box on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do you think there's a correlation? Why or why not? What other factors
%could affect glacier health and size?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
s2 = scatter(density(1:20000)*1000, depth(1:20000), 20);
xlabel('Density (kg m^{-1})')
ylabel('Depth (m)');
ylim([0 50]);
set(gca, 'Ydir', 'Reverse', 'FontSize', 14, 'LineWidth', 2);
grid on; box on;
title('Scatter Plot #2 - Global Firn Depth-Density Profiles');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do you think there's a correlation? Why or why not? What's happening in
%the shallow firn? What's happening with the high densities from 5 - 30 m
%depth?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3 Exploring the data ... beginning statistics

%% 3.1 Types of continuous probability distributions

figure(5);
subplot(1, 2, 1); %Normal
pdg = histfit(akelev, 10, 'Normal');

subplot(2, 2, 2); %Exponential
x = exprnd(700, 100, 1);
pdc = histfit(x, 10, 'exponential');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Can you think of how these distributions used for in science?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3.2 Measures of Centrality 
% Why do we use this measure? Quantifying differences in datasets can require 
% calculating some measure of central tendency. Although people commonly 
% talk about a "typical" or ?average? height, there are several standard 
% measures of ?average? (or centrality).

% The mean (also referred to as the arithmetic mean, or 
% often simply the average) is a common measure of centrality. 
% The mean is useful for symmetric distributions, but notoriously 
% sensitive to outliers. If your data set is not distributed symmetrically 
% or has extreme outliers, you will need to consider how these factors will 
% affect the calculation of the mean.
figure(5);
h = histogram(akelev, nbins);
ylabel('Number of Alaskan Glaciers');
xlabel('Mean Elevation of Accumulation Zone (m)');
set(gca, 'FontSize', 14, 'LineWidth', 2);
title('Histogram');
grid on; hold on;
plot(nanmean(akelev).*ones(1, 1800), linspace(1, 1800, 1800), 'r', ...
    'LineWidth', 2);
hold on;

% The median gives the midpoint of the sorted data, so half the data is 
% greater than the median and half is smaller. The median is much more 
% resistant than the mean to changes in a few data values, and is an 
% especially useful center for nonsymmetric (skewed) distributions, like 
% the distribution of weight data.
plot(nanmedian(akelev).*ones(1, 1800), linspace(1, 1800, 1800), 'k--', ...
    'LineWidth', 2);
hold on;
legend('mean', 'median');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In what case is the median helpful for assessing centrality?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3.3 Measures of Spread
% The difference between these extreme scenarios is the degree of spread 
% of the distributions ? that is, how much the data deviates from the 
% center. As with measures of centrality, there are several standard 
% measures of spread.

% Like the mean, the standard deviation is typically used to measure the 
% spread of symmetric distributions that follow a "bell curve" 
% (normal distribution).
% 
% Because the standard deviation is the square root of the variance ? 
% the sum of the squares of the distances of data values from the mean ? 
% the standard deviation tends to amplify the effect of outliers.
figure(6);
h = histogram(akelev, nbins);
ylabel('Number of Alaskan Glaciers');
xlabel('Mean Elevation of Accumulation Zone (m)');
set(gca, 'FontSize', 14, 'LineWidth', 2);
title('Histogram');
grid on; hold on;
plot(nanmean(akelev).*ones(1, 1800), linspace(1, 1800, 1800), 'r', ...
    'LineWidth', 2);
hold on;
plot(nanmean(akelev)-nanstd(akelev).*ones(1, 1800), ...
    linspace(1, 1800, 1800), 'r:', 'LineWidth', 2);
hold on;
plot(nanmean(akelev)+nanstd(akelev).*ones(1, 1800), ...
    linspace(1, 1800, 1800), 'r:', 'LineWidth', 2);
legend('data', 'mean', 'standard deviation');

% The interquartile range is based on the median (the 50th percentile point). 
% It gives the distance between the 25th and 75th percentile in the data ? 
% that is, the width of the region that contains the middle 50% 
% of the data values.

%% 4 Fitting Data
%% 4.1 Linear Regression 
% Suppose you suspect there is a relationship between two variables, 
% x and y. The simplest relationship (and the one you can usually assume 
% as a starting point) is that of a straight line, or y=ax+b.
% 
% The process of determining a and b for a set of x and y data is 
% called linear regression.

% The "best fit" line is the line through the data that minimizes the 
% distance between the actual, observed values of y and the values of 
% y predicted by the equation y=ax+b.
xyfit = fit(density(1:20000),depth(1:20000), 'poly1'); %Linear relationship
a = xyfit.p1;

figure(7);
s2 = scatter(density(1:20000), depth(1:20000), 20, 'k');
hold on;
p = plot(xyfit);
xlabel('Density (kg m^{-1})')
ylabel('Depth (m)');
ylim([0 65]);
set(gca, 'Ydir', 'Reverse', 'FontSize', 14, 'LineWidth', 2);
grid on; box on; hold on;
title('Linear Regression #1 - Global Firn Depth-Density Profiles');

%% 4.2 Confidence Bounds
%The fit object displays the 95% confidence bounds for each parameter. 
%These bounds mean there is a 95% likelihood that the "true" value of 
%the parameter is in that interval. In other words, if data were taken and 
%fitted repeatedly, the parameter value would be in that interval 95% of 
%the time.

%In general, the closer the 95% confidence bounds are to the value of the 
%parameter, the better the fit is.

ci = confint(xyfit, 0.95) %the 95% confidence interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What does this mean in terms of the model?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 4.3 Residuals

%Residuals are the differences between the actual data and the fit. 
%You can plot the residuals using the plot function. (The vertical distance
%from the data to the fitted line is known as the residual.)

%Residuals should be normally distributed around the zero line. If 
%there are clear outliers or a detectable pattern in the residuals, 
%the fit can be improved.
figure(8);
eqfit = 177.3*density(1:20000) + -86.66;
% yresid = depth(1:20000) - eqfit;
% SSresid = sum(yresid.^2);
% plot(density(1:20000), yresid, 'x');
plot(xyfit, density(1:20000), depth(1:20000), 'residuals');
set(gca, 'Ydir', 'Normal', 'FontSize', 14, 'LineWidth', 2);
hold on; grid on;
xlabel('Density (kg m^{-1})')
title('Residuals #1 - Global Firn Depth-Density Profiles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do you think the model fits the data well, i.e., is a linear model of 
% order 1 applicable to depth-density profiles? Do you think that firn 
% densification is a linear process?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4.4 R squared: "goodness of fit"
%Let's check your suspicions from 4.3

% R-squared is a metric of correlation. 

%Remember, ccorrelation is measured by ?r? and 
% it tells us how strongly two variables can be related. A correlation closer 
% to +1 means a strong relationship in the positive direction, while -1 means 
% a stronger relationship in the opposite direction.


% R-squared is a statistical measure of how close the data are to the fitted 
% regression line. It is also known as the coefficient of determination, or 
% the coefficient of multiple determination for multiple regression.
% 
% The definition of R-squared is fairly straight-forward; it is the 
% percentage of the response variable variation that is explained by a 
% linear model. Or:
% 
% R-squared = Explained variation / Total variation
% 
% R-squared is always between 0 and 100%:
% 
% 0% indicates that the model explains none of the variability of the 
% response data around its mean.
% 100% indicates that the model explains all the variability of the response 
% data around its mean.
% In general, the higher the R-squared, the better the model fits your data. 
% However, there are important conditions for this guideline that I?ll talk 
% about both in this post and my next post.

[xyfit, gof] = fit(density(1:20000),depth(1:20000), 'poly1');
r2 = gof.rsquare;
r2 %this is the R-squared value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How does this compare to your thoughts from section 4.3?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4.5 Linear vs. Nonlinear ... some semantics

% Linear regression can fit a nonlinear curve to the data, as long as the 
% model is linear in the coefficients, or parameters. For example, fitting a 
% cubic polynomial is still considered linear regression because the model is 
% a linear combination of the parameters p_n times a function of x.
% 
% Exponential models and sine models are two examples of nonlinear regression. 
% These models have terms that are not linear in the parameters.
% 
% For example, fitting an exponential model is nonlinear regression because 
% the parameter b is part of the exponent.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do you think our depth-density data is linear or non-linear?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's now try fitting the depth-density data with a quadratic function.

xyfit2 = fit(density(1:20000),depth(1:20000), 'poly2'); %Quadratic relationship
a = xyfit2.p1;

figure(9);
s2 = scatter(density(1:20000), depth(1:20000), 20, 'k');
hold on;
p = plot(xyfit2);
xlabel('Density (kg m^{-1})')
ylabel('Depth (m)');
ylim([0 65]);
set(gca, 'Ydir', 'Reverse', 'FontSize', 14, 'LineWidth', 2);
grid on; box on; hold on;
title('Linear Regression #2 - Global Firn Depth-Density Profiles');

%the fit can be improved.
figure(10);
eqfit = 273.4.*density(1:20000).^2 + -172.*density(1:20000) + 21.94; %change this if you 
%change the data or fit
plot(xyfit2, density(1:20000), depth(1:20000), 'residuals');
set(gca, 'Ydir', 'Normal', 'FontSize', 14, 'LineWidth', 2);
hold on; grid on;
xlabel('Density (kg m^{-1})')
title('Residuals #2 - Global Firn Depth-Density Profiles');

[xyfit, gof] = fit(density(1:20000),depth(1:20000), 'poly2');
r2_2 = gof.rsquare;
r2_2 %this is the R-squared value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the fit improve? How do you know?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 4.6 A model
% How do we use measurements to inform our decisions about best-fit 
% equations to certain processes? Let's look at an empirical model for 
% firn densification that uses density data. It uses temperature and
% accumulation rate as the "boundary conditions" and empirically tuned
% constants to fit a curve to the depth-density data.

%The empiricaly based model:
T = (273+(-20)); % firn temperature (Kelvin) 
acc = 0.45;% ICE-equiv accumulation rate (m/year) 
SEP = 0;% distance Tx to Rx (metres)   
c = 3e8; c_ice = 1.685e8; E_ice = (c/c_ice)^2;  
rho_i = 0.917; 
rho_0 = 0.4;% firn density at surface 
rho_c = 0.55; % critical density (stage 1-2 densification transition) 
rho_w = 1; % water density    
% density  vs. depth profile (look-up table to 5000 metres) 
% from Herron-Langway model 
% rho vs. z, where rho is in Mg per cubic metres 
z = [0:0.05:5000]'; 
R = 8.314; 
k0 = 11*exp(-10160/(R*T)); 
k1 = 575*exp(-21400/(R*T)); 
alpha = log(rho_0/(rho_i-rho_0)); 
beta = log(rho_c/(rho_i-rho_c));  
A = acc*rho_i/rho_w;% WATER-equiv accumulation rate (m/year) 
hc = (beta-alpha)/(rho_i*k0); % critical depth hc (metre)   
Z=(z<hc).*exp(rho_i*k0*z+alpha) + ... % density profile evaluation
    (z>=hc).*exp(rho_i*k1*(z-hc)/(A^0.5)+beta); 
rho=rho_i*Z./(1+Z);   

figure(11);
s2 = scatter(density(1:20000), depth(1:20000), 20, 'k');
hold on;
p = plot(rho, z, 'r', 'LineWidth', 2);
xlabel('Density (kg m^{-1})')
ylabel('Depth (m)');
ylim([0 65]);
set(gca, 'Ydir', 'Reverse', 'FontSize', 14, 'LineWidth', 2);
grid on; box on; hold on;
title('Empirically Based Model and Global Firn Depth-Density Profiles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What do you think about this model? Do you think that this model was 
%tuned to depth-density data across many temperature and accumulation 
% conditions? (Hint: we had to use an accumulation rate of 4.5 m per year
% as the boundary condition ... do you think this is representative of
% where the global firn cores were retrieved?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

