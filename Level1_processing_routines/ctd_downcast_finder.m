function [ind_start_down, ind_end_down, ind_start_up, ind_end_up, counter, dpdt_inst] = ctd_downcast_finder(P, time, castnum)
%% This function identifies downcasts from the EcoCTD data

% This function is called by ctd_lag_finder.m, plot_eco_casts.m, and rsk2netcdf.m.
% This function finds downcasts and upcasts in EcoCTD file identified by
% castnum (the serial number for the cast).

% compute fall rate (using first differencing)
dpdt = cat(1,0,diff(P)./diff(time)/86400);
% compute a smoothed version
dpdtlow = smooth(dpdt,80);
% Set a velocity threshold to define free-falling
velocity_threshold = 1;

% Calculate more precise version of dpdt using center differencing.
% First smooth pressure using low-pass filter, cutoff at 2 sec
[P_smooth,~,~] = bwfilter(P,2/3600,5,8);
dpdt_inst = NaN(length(P_smooth),1);
for kk = 1:length(P_smooth)
    % decide on interval width
    h = 2/8; % use 2 scans (convert to seconds)
    % second-order forward difference
    if kk == 1
        dpdt_inst(kk) = (4*P_smooth(2)-P_smooth(3)-3*P_smooth(1))/(2*h);
    % second-order backward difference
    elseif kk == length(P_smooth)
        dpdt_inst(kk) = (-4*P_smooth(end-1)+1*P_smooth(end-2)+3*P_smooth(end))/(2*h);
    % first-order center difference
    else
        dpdt_inst(kk) = (P_smooth(kk+1) - P_smooth(kk-1))/h;
    end
end

% Start a counter for the number of profiles within the file
% Set initial start point to look for downcasts at 1
counter = 0;
start_pt = 1;
% loop through P values, looking for start and end points of downcasts and upcasts
% continue loop as long as we can find a start point
% start points are defined by points close to the surface (between 1.5m
% and 5m) where the speed of the probe crosses the velocity threshold.
while (isempty(find(P(start_pt:end)>1.5 & P(start_pt:end)<5 ...
        & dpdtlow(start_pt:end)>velocity_threshold,1,'first')))==0
    % initialize flags for extra profiles or bad profiles
    extra_profile = 0;
    bad_profile = 0;
    % find starting point of downcast
    ind_start_down(counter+1) = start_pt + ...
        find(P(start_pt:end)>1.5 & P(start_pt:end)<5 ...
        & dpdtlow(start_pt:end)>velocity_threshold,1,'first');
    % check to make sure start point isn't the last point of the file
    if ind_start_down(counter+1) == length(P)
        ind_start_down(counter+1) = [];
        break
    end
    % find the end point -- where depth is deeper than 1.5m and the
    % (smoothed) velocity of the probe changes from positive 
    % (going down) to negative (coming back up)
    ind_end_down(counter+1) = ind_start_down(counter+1) + ...
        find(P(ind_start_down(counter+1):end)>1.5 & ...
        dpdtlow(ind_start_down(counter+1):end)<=0,1,'first') - 1;
    % if we overshot the end of the profile, track back by finding
    % where P stopped increasing
    if P(ind_end_down(counter+1))-P(ind_end_down(counter+1)-1) < 0
        ind_end_down(counter+1) = ind_end_down(counter+1) - find(flipud(diff(P(ind_start_down(counter+1):ind_end_down(counter+1))))>0,1,'first')+1;
    end
    
    % get upcast
    ind_start_up(counter+1) = ind_end_down(counter+1)+1;
    ind_end_up(counter+1) = ind_start_up(counter+1) + ...
        find(P(ind_start_up(counter+1):end)<1.5 & ...
        dpdtlow(ind_start_up(counter+1):end)>=-1*velocity_threshold,1,'first') - 1;                

    % check for extra short profiles (which we throw out) or profiles 
    % that contain a segment of decreasing depth somewhere in the middle
    % these should be split into two pieces to get rid of extra points
    % when the probe was moving upwards
    if P(ind_end_down(counter+1)) < 200
        % if the profile ends before max depth, check to see if it is
        % closely followed by another profile start point, within the
        % next 10 seconds
        blip = find(dpdtlow(ind_end_down(counter+1)+8:ind_end_down(counter+1)+88)>velocity_threshold,1,'first');
        if blip
            % if next profile starts near the surface, write over the
            % previous start and end points using new complete profile
            % (this implies the probe went down a bit, then came back up
            % and started a new profile -- so we can throw out the blip)
            if P(ind_end_down(counter+1)+blip) > 1.5 && P(ind_end_down(counter+1)+blip) < 5
                ind_start_down(counter+1) = ind_end_down(counter+1) + blip;
                ind_end_down(counter+1) = ind_start_down(counter+1) + ...
                    find(P(ind_start_down(counter+1):end)>1.5 & ...
                    dpdtlow(ind_start_down(counter+1):end)<=0,1,'first') - 1;
                if P(ind_end_down(counter+1))-P(ind_end_down(counter+1)-1) < 0
                    ind_end_down(counter+1) = ind_end_down(counter+1) - find(flipud(diff(P(ind_start_down(counter+1):ind_end_down(counter+1))))>0,1,'first')+1;
                end
                
                % get upcast
                ind_start_up(counter+1) = ind_end_down(counter+1)+1;
                ind_end_up(counter+1) = ind_start_up(counter+1) + ...
                    find(P(ind_start_up(counter+1):end)<1.5 & ...
                    dpdtlow(ind_start_up(counter+1):end)>-1*velocity_threshold,1,'first') - 1;
            % if next profile doesn't start near the surface, add a new
            % profile to get a more complete picture (keep both sections)
            else
                ind_start_down(counter+2) = ind_end_down(counter+1) + blip;
                ind_end_down(counter+2) = ind_start_down(counter+2) ...
                    + find(P(ind_start_down(counter+2):end)>1.5 & ...
                    dpdtlow(ind_start_down(counter+2):end)<=0,1,'first') - 1;
                extra_profile = 1;
                if P(ind_end_down(counter+2))-P(ind_end_down(counter+2)-1) < 0
                    ind_end_down(counter+2) = ind_end_down(counter+2) - find(flipud(diff(P(ind_start_down(counter+2):ind_end_down(counter+2))))>0,1,'first')+1;
                end
                
                % get both upcasts
                ind_start_up(counter+1) = ind_end_down(counter+1)+1;
                ind_end_up(counter+1) = ind_start_up(counter+1) + ...
                    find(dpdtlow(ind_start_up(counter+1):end)>-1*velocity_threshold,1,'first') - 1;
                ind_start_up(counter+2) = ind_end_down(counter+2)+1;
                ind_end_up(counter+2) = ind_start_up(counter+2) + ...
                    find(P(ind_start_up(counter+2):end)<1.5 & ...
                    dpdtlow(ind_start_up(counter+2):end)>-1*velocity_threshold,1,'first') - 1;
            
            end
        % if not closely followed by another profile
        else
            % if profile is too short to count on its own, throw it out
            if P(ind_end_down(counter + 1)) - P(ind_start_down(counter+1)) < 10
                % save point after end of bad profile so we can skip it
                % during the next iteration of the while loop
                new_start = ind_end_down(counter+1);
                ind_end_down(counter+1) = [];
                ind_start_down(counter+1) = [];
                bad_profile = 1;
                % no upcast for bad profile
            % otherwise it's just a normal profile, count as usual
            % note: there will be issues if the find function returns empty
            % vectors, because this code tries to add find to ind_start!
                
            end
        end
            
    end
    % if we didn't find an end point, (ind_end is empty), use last point in file
    if bad_profile == 0 && isempty(ind_end_up)
        ind_end_up = length(P);
    end
    % use flags to determine appropriate counter and start point for
    % next iteration of loop
    if extra_profile
        counter = counter+2;
        start_pt = ind_end_up(counter);
    elseif bad_profile
        start_pt = new_start;
    else
        counter = counter + 1;
        start_pt = ind_end_up(counter);
    end
end

disp(['This is file number ',num2str(castnum)])
disp(['This file contains ',num2str(counter),' profiles'])
end
