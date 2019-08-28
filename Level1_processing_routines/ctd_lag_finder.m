% This code gets all the upcasts and downcasts from the EcoCTD data and 
% calculates the lag correction to minimize salinity spiking. 
% We calculate the lag that maximizes the correlation between 
% first differences in T and C. We do this for 60 scan segments of every 
% profile, then determine the mode of the resulting distribution.
clear
% files we want to look at
filepath = './'
files = dir([filepath 'Level0/*.rsk']);

P_num = 1;
seg_num = 1;
% Save lag for each segment in a matrix (rows = number of profiles, 
% columns = number of 60-scan segments in each profile
% use preallocated matrix of NaNs to avoid padding with extra zeros, which alters distribution of lags
lag = NaN(113,73); xc = NaN(113,73); %fall_rate = NaN(1697,38);
lag_up = NaN(113,73); xc_up = NaN(113,73);
lag_down = NaN(113,73); xc_down = NaN(113,73);

% Loop through each RSK file in the folder.
for FF = 1:length(files)
        
    % Generate filename
    RSKfilename = fullfile(files(FF).folder,files(FF).name);
    
    % Open the RSK file using RBR Ruskin toolbox
    RSK = RSKopen(RSKfilename);
    % Extract the data from RSK file
    RSK = RSKreaddata(RSK);

    % Record Section number from filename
    disp(['Section number: ',RSKfilename(end-5:end-4)]);
    Sec_num = str2double(RSKfilename(end-5:end-4));

    % Get time, pressure, conductivity, temperature
    time = RSK.data.tstamp(:);
    P = RSK.data.values(:,3)-10.1325; %correct for atmospheric pressure
    T = RSK.data.values(:,2);
    C = RSK.data.values(:,1);
    
    % Find downcasts and number of profiles (one down and one up = one profile)
    [ind_start_down, ind_end_down, ind_start_up, ind_end_up,counter, dpdt] = ctd_downcast_finder(P, time, Sec_num);


    % Calculate maximum cross-covariance of first differences in T and C
    % and store lag values in a list.
    % Before calculating lag, split profiles into 60 scan segments.
    % Leftover scans at the bottom of each profile are ignored.
    % Fit a 2nd order polynomial to the lags for each segment and save the
    % lag of maximum correlation for each profile.
    % Then find the mode of the resulting distribution of lags (the lag of
    % maximum correlation for the entire distribution). 
    
    % Loop through profiles in this data file
    for ii = 1:length(ind_start_down)   
        % Find data indices for the considered profile, for both upcasts
        % and downcasts.
        %ind_down = find(RSK.data.tstamp>=RSK.profiles.downcast.tstart(ii) &...
        %    RSK.data.tstamp<=RSK.profiles.downcast.tend(ii));
        %ind_up = find(RSK.data.tstamp>=RSK.profiles.upcast.tstart(ii) &...
        %    RSK.data.tstamp<=RSK.profiles.upcast.tend(ii)); 
        ind_down = ind_start_down(ii):ind_end_down(ii);
        ind_up = ind_start_up(ii):ind_end_up(ii);
        
        % Analyze both downcast and upcast
        direction = {'down','up'};
%        seg_num = 1;
        for dd = 1:length(direction)
            if strcmp(direction{dd},'up')
                theind = ind_up;
            elseif strcmp(direction{dd},'down')
                theind = ind_down;
            end
            % get temperature and conductivity data
            temp = T(theind);
            cond = C(theind);
            seg_num = 1;
            % split each profile into 60-scan segments
            % record lag of max covariance and average dpdt over this segment
            for pp = 1:floor(length(temp)/60)
                seg_index = (pp-1)*60+1:pp*60; 
                % get cross-covariance of first differences and store
                [c,lags] = xcov(diff(temp(seg_index)), diff(cond(seg_index)));
                % find lag with maximum covariance
                [~,I] = max(abs(c));
                % fit a polynomial to 3 points around the max for this segment
                % (exceptions if the max happens to be at an edge)
                if I == length(lags)
                    p = polyfit(lags(I-2:I)',c(I-2:I),2);
                elseif I == 1    
                    p = polyfit(lags(I:I+2)',c(I:I+2),2);
                else   
                    p = polyfit(lags(I-1:I+1)',c(I-1:I+1),2);
                end
                % get the coefficients of the polynomial
                a = p(1);
                b = p(2);
                % calculate and save the lag and xcov where the maximum occurs    
                max_lag = -b/(2*a);
                if strcmp(direction{dd},'up')
                    lag_up(P_num,seg_num) = max_lag;
                    xc_up(P_num,seg_num) = polyval(p,max_lag);
                elseif strcmp(direction{dd},'down')
                    lag_down(P_num,seg_num) = max_lag;
                    xc_down(P_num,seg_num) = polyval(p,max_lag);
                end
                % update segment number
                seg_num = seg_num + 1;    
            end        
        end         
        % update profile number
        P_num = P_num + 1;
    end; clear ii
end; clear FF

%% Plot distribution of lags
figure(1)
clf
% first plot entire distribution
%subplot(2,1,1)
% histogram(lag(:),-10:.2:10)
% xlabel('lag of maximum cross covariance')
% ylabel('count')
% title('full distribution')
% then do an iterative cut to get rid of tails
%{
mean_lag = nanmean(lag(:));
stdev_lag = nanstd(lag(:));
center_ind = find((lag > mean_lag-1*stdev_lag) & (lag < mean_lag+1*stdev_lag));
new_lag = lag(center_ind);
new_fall = fall_rate(center_ind);
for ii = 1:3
    mean_lag = nanmean(new_lag);
    stdev_lag = nanstd(new_lag);
    center_ind = find((new_lag > mean_lag-2*stdev_lag) & (new_lag < mean_lag+2*stdev_lag));
    new_lag = new_lag(center_ind);
    new_fall = new_fall(center_ind);
end

subplot(2,1,2)
histogram(new_lag,-10:.2:10)
xlabel('lag of maximum cross covariance')
ylabel('count')
title('distribution after cut')
%}
% make a fit for the center of the distribution and calculate max
% [N,edges] = histcounts(lag, -2:.1:5);
% [~,I] = max(N);
% max_lag = (edges(I)+edges(I+1))/2;
% disp(['The lag of maximum correlation is: ', num2str(max_lag), ' scans.'])

subplot(2,1,1)
histogram(lag_up(:),-2:.05:3)
xlabel('lag of maximum cross covariance')
ylabel('count')
title('upcasts')
subplot(2,1,2)
histogram(lag_down(:),-2:.05:3)
xlabel('lag of maximum cross covariance')
ylabel('count')
title('downcasts')

%% Make 3D histogram showing lags and covariances
figure(2)
subplot(2,1,1)
hist3([lag_down(:), xc_down(:)],'Edges',{-2:.05:3 -.08:.001:.08})
xlabel('lag of maximum covariance (number of scans)')
ylabel('maximum covariance')
zlabel('count')
title('EcoCTD lag covariances, downcast')

subplot(2,1,2)
hist3([lag_up(:), xc_up(:)],'Edges',{-2:.05:3 -.08:.001:.08})
xlabel('lag of maximum covariance (number of scans)')
ylabel('maximum covariance')
zlabel('count')
title('EcoCTD lag covariances, upcast')

%% Calculate a weighted average to find max lag
% use covariance to weight lag values
wei_mean_down = sum(lag_down(:).*xc_down(:),'omitnan')/sum(xc_down(:),'omitnan');
mean_down = nanmean(lag_down(:));
median_down = median(lag_down(:),'omitnan');

%% Plot distribution of lags vs. dpdt
%{
figure(2)
clf
plot(fall_rate(:), lag(:), 'k.', 'MarkerSize',2)
hold on
plot(new_fall, new_lag, 'b.','MarkerSize',5)
xlabel('fall rate (dbar/s)')
ylabel('lag (number of scans)');
% fit a linear regression to look for a trend

% need to choose data range before performing regression (lots of outliers)
% create bins for lag values in 0.25 dbar/s bins (as in Ullman and Hebert 2014)
lag_bins = 1.375:.25:3.625;
good_lag = [];
good_fall = [];
mean_lag = NaN(1,length(lag_bins)-1);
mean_fall = NaN(1,length(lag_bins)-1);
% find average lag in each bin
for ii = 1:(length(lag_bins) - 1)
    bin_start = lag_bins(ii);
    bin_end = lag_bins(ii+1);
    bin_ind = find(new_fall>=bin_start & new_fall<bin_end);
    binned_fall = new_fall(bin_ind);
    binned_lag = new_lag(bin_ind);
    bin_mean = nanmean(binned_lag);
    % find standard deviation, remove points greater than 2sd from mean, recalculate mean
    bin_stdev = nanstd(binned_lag);
    new_ind = find((binned_lag >= bin_mean-1*bin_stdev) & (binned_lag <= bin_mean+1*bin_stdev));
    % save mean lag and all points within the criteria
    good_lag = [good_lag; binned_lag(new_ind)];
    mean_lag(ii) = mean(binned_lag(new_ind));
    mean_fall(ii) = (bin_start+bin_end)/2;
    good_fall = [good_fall; binned_fall(new_ind)];
end; clear ii

% Plot points that met criteria plus means
figure(3)
clf
plot(good_fall, good_lag,'b.','MarkerSize',5)
hold on
plot(mean_fall, mean_lag, 'r.-','MarkerSize',10)
xlabel('fall rate (dbar/s)')
ylabel('lag (number of scans)');

%% Get an equation for that curve so we can get the lag for any given fall rate
% Use model developed by Ullman and Hebert 2014

% first define constants used in model (physical parameters describing probe)
L_cell = 0.11; % m
delta_x = 0.021; % m
l = 0.17; % m
v = 1.36e-6; % m^2/s
a = 2e-3; % m

% get values for u2 from u1 = dpdt = bins calculated above
u1 = mean_fall;
u2 = (-8*v*l + sqrt((8*v*l)^2 + a^4*u1.^2))/(a^2);
tau = mean_lag;

% then get a fit for our binned average curve -- using the equation
% tau = c0 + c2/u2 
y = tau';
x = u2';
ft = fittype({'1','1./x'});

lag_model = fit(x,y,ft);

% get the constants we need
c0 = lag_model.a;
c2 = lag_model.b;

% plot
figure(4)
clf
plot(lag_model,x,y)
%}
