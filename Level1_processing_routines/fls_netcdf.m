%% This code generates the Level-1 NetCDF files for the EcoCTD.
% It creates a file for downcasts and upcasts for the variables
% measured -- and computed from -- the ECOPuck.
%
%
filepath = './'
%% Processing steps are
%       1- Glitch removal.
%           This was attributed to "Analog to Digital (A2D)
%           calibration, which repeats previously logged values. Exact
%           repetition of a datapoint in pressure, chl, or backscatter
%           are removed from the data.
%
%       2- Apply lag correction
%           The lag correction to pair CTD measurements with ECOPuck takes
%           2 forms: (1) the lag resulting from the physical location of
%           the sensors on the instrument (72 cm apart on the EcoCTD), and
%           the lag due to the response time of the intrument. The former
%           is fall-rate-dependent, while the other is constant and defined
%            [in number of scans]
lag = 0;
%
%       3- Correction for atmospheric pressure (-10.1325)
%
%       4- Apply cross-calibration factors (y = SF*x + offset) determined
%       from ship-based CTD.
SF = 0.00446;
offset = 0.44634;
%
%       5- Set the values for flag

direction = {'down','up'};

for dd = 1:length(direction)
    
    % Select the proper indices based on cast direction (up or down)
    if strcmp(direction{dd},'up')
        theind = ind_up;
    elseif strcmp(direction{dd},'down')
        theind = ind_down;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Process ECOPuck data
    %%%%%%%%%%%%%%%%%%%%%
    
    %% Step 0 - extract variables from RSK file
    % Pressure
    FLS_pres = RSK.data.values(theind,3);
    % Time
    % Correct for date issue on Section 16
    if Sec_num == 16
        FLS_time = RSK.data.tstamp(theind)+6724.46374999;
    else
        FLS_time = RSK.data.tstamp(theind);
    end
    FLS_bb470 = RSK.data.values(theind,5);
    FLS_bb700 = RSK.data.values(theind,6);
    FLS_chl = RSK.data.values(theind,7);
    
    %% STEP 1 - Glitch removal. aka "Zero-order hold correction"
    % Must test different channels separately.
    %
    %     dpdt = diff(FLS_pres); rmv = find(dpdt == 0)+1; FLS_pres(rmv) = NaN;
    %     dtdt = diff(FLS_time); rmv = find(dtdt == 0)+1; FLS_time(rmv) = NaN;
    %     dbb470dt = diff(FLS_bb470); rmv = find(dbb470dt == 0)+1; FLS_bb470(rmv) = NaN;
    %     dbb700dt = diff(FLS_bb700); rmv = find(dbb700dt == 0)+1; FLS_bb700(rmv) = NaN;
    %     dchldt = diff(FLS_chl); rmv = find(dchldt == 0)+1; FLS_chl(rmv) = NaN;
%     dpdt = diff(FLS_pres); dbb470dt = diff(FLS_bb470); dbb700dt = diff(FLS_bb700); dchldt = diff(FLS_chl);
%     ind = find(dpdt == 0 | dbb470dt == 0 | dbb700dt == 0 | dchldt == 0);
%     FLS_pres(ind+1) = NaN;
%     FLS_bb470(ind+1) = NaN;
%     FLS_bb700(ind+1) = NaN;
%     FLS_chl(ind+1) = NaN;
%     clear dpdt dtdt dbb470dt dbb700dt dchldt rmv ind
%     
    %% STEP 2 - Apply lag correction
    
    % First, compute the lag due to physical separation of the sensors
    % Fall rate
    ws = cat(1,0,diff(FLS_pres)./(diff(FLS_time)*86400));
    % Lag in number of scans (sensor separation = 0.72 m, 8Hz sampling)
    lagadv = 0.72./ws*8;
    
    % Second, apply the lags
    FLS_bb470(:,1) = interp1(1:length(FLS_bb470),FLS_bb470,(1:length(FLS_bb470))+lagadv');
    FLS_bb700(:,1) = interp1(1:length(FLS_bb700),FLS_bb700,(1:length(FLS_bb700))+lagadv');
    FLS_chl(:,1) = interp1(1:length(FLS_chl),FLS_chl,(1:length(FLS_chl))+lagadv');
    
    %% STEP 3 - Correction for atmospheric pressure
    FLS_SeaPress = FLS_pres-10.1325;
    
    %% STEP 4 - Secondary variables
    % Longitude & Latitude
    % Get boat position at each timestamp
    [lon,lat] = findcoordinates(FLS_time);
    % Only uses first value (i.e. dropping location) for (lon,lat)
    if strcmp(direction{dd},'down')
        theFLSlondown = lon(1);
        theFLSlatdown = lat(1);
        FLS_lon = lon(1)*ones(length(theind),1); clear lon
        FLS_lat = lat(1)*ones(length(theind),1); clear lat
        % Linear interpolation in time between drop point and boat location
        % when put on board.
    elseif strcmp(direction{dd},'up')
        FLS_lon = interp1([FLS_time(1) FLS_time(end)],[theFLSlondown lon(end)],FLS_time);
        FLS_lat = interp1([FLS_time(1) FLS_time(end)],[theFLSlatdown lat(end)],FLS_time);
        clear theFLSlondown theFLSlatdown
    end
    
    %% Step 5 - Apply cross-calibration factors
    
    FLS_chlugL = SF*FLS_chl+offset;
    
    %% Step 6 - Set the value flag
    FLS_bb470_flag = ones(length(FLS_SeaPress),1);
    FLS_bb700_flag = ones(length(FLS_SeaPress),1);
    FLS_chlugL_flag = ones(length(FLS_SeaPress),1);
    
    % Put a zero Flag is drop rate is too slow
    FLS_bb470_flag(abs(ws)<1) = 0;
    FLS_bb700_flag(abs(ws)<1) = 0;
    FLS_chlugL_flag(abs(ws)<1) = 0;
    
    % Zero flag if less than 10 m on the way down, 30 m on the way up
    if strcmp(direction,'up')
        FLS_bb470_flag(FLS_SeaPress<30) = 0;
        FLS_bb700_flag(FLS_SeaPress<30) = 0;
        FLS_chlugL_flag(FLS_SeaPress<30) = 0;
    elseif strcmp(direction,'down')
        FLS_bb470_flag(FLS_SeaPress<10) = 0;
        FLS_bb700_flag(FLS_SeaPress<10) = 0;
        FLS_chlugL_flag(FLS_SeaPress<10) = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Generate NetCDF file
    %%%%%%%%%%%%%%%%%%%%%
    
    % Generate filename
    filename = [filepath 'EcoCTD_fls/EcoCTD_fls_S',num2str(Sec_num,'%02d'),...
        '_P',num2str(P_num,'%03d'),'_',direction{dd},'.nc'];
    
    % Create the NetCDF file
    ncid = netcdf.create(filename,'NOCLOBBER');
    
    % Define dimensions
    time_dimID = netcdf.defDim(ncid,'Time',length(theind));
    
    % Define Global Attributes
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'title','Level-1 processed data from EcoPuck instrument on EcoCTD package');
    netcdf.putAtt(ncid,varid,'institution','Woods Hole Oceanographic Institution');
    netcdf.putAtt(ncid,varid,'source','ocean profile observations');
    netcdf.putAtt(ncid,varid,'history',[datestr(now),' - File generated by Dr. M. Dever']);
    netcdf.putAtt(ncid,varid,'references','CALYPSO Data Report');
    netcdf.putAtt(ncid,varid,'external_variables','');
    netcdf.putAtt(ncid,varid,'Conventions','CF-1.7');
    netcdf.putAtt(ncid,varid,'creation_date',datestr(now));
    netcdf.putAtt(ncid,varid,'Comment','');
    
    % Define variables and their attributes
    Pid = netcdf.defVar(ncid,'P','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,Pid,'long_name','Absolute Pressure');
    netcdf.putAtt(ncid,Pid,'standard_name','sea_water_pressure');
    netcdf.putAtt(ncid,Pid,'units','dbar');
    netcdf.putAtt(ncid,Pid,'valid_range',[0 600]);
    netcdf.putAtt(ncid,Pid,'actual_range',[min(FLS_pres(:)) max(FLS_pres(:))]);
    netcdf.putAtt(ncid,Pid,'missing_value',-999);
    netcdf.putAtt(ncid,Pid,'origin','Measured');
    
    tid = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,tid,'long_name','Time');
    netcdf.putAtt(ncid,tid,'units','days');
    netcdf.putAtt(ncid,tid,'axis','T');
    netcdf.putAtt(ncid,tid,'valid_range',[min(FLS_time(:)) max(FLS_time(:))]);
    netcdf.putAtt(ncid,tid,'actual_range',[min(FLS_time(:)) max(FLS_time(:))]);
    netcdf.putAtt(ncid,tid,'missing_value',-999);
    netcdf.putAtt(ncid,tid,'origin','Measured');
    netcdf.putAtt(ncid,tid,'notes','days since [0000 0 0 0 0 0]');
    
    lonid = netcdf.defVar(ncid,'lon','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,lonid,'long_name','Longitude');
    netcdf.putAtt(ncid,lonid,'standard_name','longitude');
    netcdf.putAtt(ncid,lonid,'units','degree_east');
    netcdf.putAtt(ncid,lonid,'valid_range',[-180 180]);
    netcdf.putAtt(ncid,lonid,'actual_range',[min(FLS_lon(:)) max(FLS_lon(:))]);
    netcdf.putAtt(ncid,lonid,'missing_value',-999);
    netcdf.putAtt(ncid,lonid,'origin','Computed');
    
    latid = netcdf.defVar(ncid,'lat','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,latid,'long_name','Latitude');
    netcdf.putAtt(ncid,latid,'standard_name','latitude');
    netcdf.putAtt(ncid,latid,'units','degree_north');
    netcdf.putAtt(ncid,latid,'valid_range',[-90 90]);
    netcdf.putAtt(ncid,latid,'actual_range',[min(FLS_lat(:)) max(FLS_lat(:))]);
    netcdf.putAtt(ncid,latid,'missing_value',-999);
    netcdf.putAtt(ncid,latid,'origin','Computed');
    
    sPid = netcdf.defVar(ncid,'seaPress','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,sPid,'long_name','Sea Pressure');
    netcdf.putAtt(ncid,sPid,'standard_name','sea_water_pressure_due_to_sea_water');
    netcdf.putAtt(ncid,sPid,'units','dbar');
    netcdf.putAtt(ncid,sPid,'valid_range',[0 600]);
    netcdf.putAtt(ncid,sPid,'actual_range',[min(FLS_SeaPress(:)) max(FLS_SeaPress(:))]);
    netcdf.putAtt(ncid,sPid,'missing_value',-999);
    netcdf.putAtt(ncid,sPid,'origin','Computed');
    netcdf.putAtt(ncid,sPid,'notes','Computed as P-10.1325 dbar');
    
    bb470id = netcdf.defVar(ncid,'bb470','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,bb470id,'long_name','Backscatter at 470 nm');
    netcdf.putAtt(ncid,bb470id,'units','1');
    netcdf.putAtt(ncid,bb470id,'valid_range',[0 max(FLS_bb470(:))]);
    netcdf.putAtt(ncid,bb470id,'actual_range',[min(FLS_bb470(:)) max(FLS_bb470(:))]);
    netcdf.putAtt(ncid,bb470id,'missing_value',-999);
    netcdf.putAtt(ncid,bb470id,'ancillary_variables','bb470_qc');
    netcdf.putAtt(ncid,bb470id,'origin','Measured');
    netcdf.putAtt(ncid,bb470id,'notes','Non-standard units of [counts]');
    
    bb470flagid = netcdf.defVar(ncid,'bb470_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,bb470flagid,'long_name','Backscatter at 470 nm Quality');
    netcdf.putAtt(ncid,bb470flagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,bb470flagid,'units','1');
    netcdf.putAtt(ncid,bb470flagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,bb470flagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,bb470flagid,'missing_value',-999);
    netcdf.putAtt(ncid,bb470flagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,bb470flagid,'flag_meanings','bad questionable good');
    
    bb700id = netcdf.defVar(ncid,'bb700','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,bb700id,'long_name','Backscatter at 700 nm');
    netcdf.putAtt(ncid,bb700id,'units','1');
    netcdf.putAtt(ncid,bb700id,'valid_range',[0 max(FLS_bb700(:))]);
    netcdf.putAtt(ncid,bb700id,'actual_range',[min(FLS_bb700(:)) max(FLS_bb700(:))]);
    netcdf.putAtt(ncid,bb700id,'missing_value',-999);
    netcdf.putAtt(ncid,bb700id,'ancillary_variables','bb700_qc');
    netcdf.putAtt(ncid,bb700id,'origin','Measured');
    netcdf.putAtt(ncid,bb700id,'notes','Non-standard units of [counts]');
    
    bb700flagid = netcdf.defVar(ncid,'bb700_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,bb700flagid,'long_name','Backscatter at 700 nm Quality');
    netcdf.putAtt(ncid,bb700flagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,bb700flagid,'units','1');
    netcdf.putAtt(ncid,bb700flagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,bb700flagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,bb700flagid,'missing_value',-999);
    netcdf.putAtt(ncid,bb700flagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,bb700flagid,'flag_meanings','bad questionable good');
    
    chlid = netcdf.defVar(ncid,'chl_raw','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,chlid,'long_name','Raw Chlorophyll-Fluorescence');
    netcdf.putAtt(ncid,chlid,'units','1');
    netcdf.putAtt(ncid,chlid,'origin','Measured');
    
    chlugLid = netcdf.defVar(ncid,'chl_cal','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,chlugLid,'long_name','Calibrated Chlorophyll-Fluorescence');
    netcdf.putAtt(ncid,chlugLid,'units','microgram liter-1');
    netcdf.putAtt(ncid,chlugLid,'valid_range',[0 max(FLS_chl(:))]);
    netcdf.putAtt(ncid,chlugLid,'actual_range',[min(FLS_chl(:)) max(FLS_chl(:))]);
    netcdf.putAtt(ncid,chlugLid,'missing_value',-999);
    netcdf.putAtt(ncid,chlugLid,'ancillary_variables','chl_cal_qc');
    netcdf.putAtt(ncid,chlugLid,'origin','Computed');
    
    chlflagid = netcdf.defVar(ncid,'chl_cal_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,chlflagid,'long_name','Calibrated Chlorophyll-Fluorescence Quality');
    netcdf.putAtt(ncid,chlflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,chlflagid,'units','1');
    netcdf.putAtt(ncid,chlflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,chlflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,chlflagid,'missing_value',-999);
    netcdf.putAtt(ncid,chlflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,chlflagid,'flag_meanings','bad questionable good');
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Write Data to Netcdf File
    %%%%%%%%%%%%%%%%%%%%%
    
    % Absolute Pressure
    netcdf.putVar(ncid,Pid,FLS_pres);
    
    % Time
    netcdf.putVar(ncid,tid,FLS_time);
    
    % Sea Pressure
    netcdf.putVar(ncid,sPid,FLS_SeaPress);
    
    % backscatter (470nm)
    netcdf.putVar(ncid,bb470id,FLS_bb470);
    
    % backscatter (470nm) flag
    netcdf.putVar(ncid,bb470flagid,FLS_bb470_flag);
    
    % backscatter (700nm)
    netcdf.putVar(ncid,bb700id,FLS_bb700);
    
    % backscatter (700nm) flag
    netcdf.putVar(ncid,bb700flagid,FLS_bb700_flag);
    
    % Chlorophyll raw
    netcdf.putVar(ncid,chlid,FLS_chl);
    
    % Chlorophyll calibrated
    netcdf.putVar(ncid,chlugLid,FLS_chlugL);
    
    % Chlorophyll calibrated
    netcdf.putVar(ncid,chlflagid,FLS_chlugL_flag);
    
    % Longitude & Latitude
    netcdf.putVar(ncid,lonid,FLS_lon);
    netcdf.putVar(ncid,latid,FLS_lat);
    
    netcdf.close(ncid);
    clear filename theind ncid *ID *id FLS*
    
end
