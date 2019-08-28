%% This code generates the Level-1 NetCDF files for the EcoCTD.
% It creates a file for downcasts and upcasts for the variables
% measured -- and computed from -- the oxygen probe.
%
%
filepath = './'
%% Processing steps are
%       1- Glitch removal.
%           This was attributed to "Analog to Digital (A2D)
%           calibration, which repeats previously logged values. Exact
%           repetition of a datapoint in pressure, time, or oxygen
%           saturation are removed from the data.
%
%       2- Apply lag correction
%           The lag correction to pair CTD measurements with optode takes
%           2 forms: (1) the lag resulting from the physical location of
%           the sensors on the instrument (48 cm apart on the EcoCTD), and
%           the lag due to the response time of the intrument. The former
%           is fall-rate-dependent, and determined in the core code, while
%           the other is constant and defined below [in number of scans].
lag = 6;
%
%       3- Correction for atmospheric pressure (-10.1325)
%
%       4- Secondary variables are computed using the TEOS-10 Gibbs
%           Seawater package. These variables are:
%               - Oxygen concentration (in mL/L)
%               - Oxygen concentration (in umol/L)
%               - Longitude and Latitude. On downcast: boat position at
%               drop point. On upcast: linear interpolation in time between
%               drop location at beginning of upcast, and boat location at
%               end of upcast.
%
%       5- Apply cross-calibration factors (y = SF*x + offset)
SF = 0.94356;offset = 0.00057;
%
%       6 - Applies a 16 point median filter to oxygen (~2s)
%
%       7- Set the values for flag

direction = {'down','up'};

for dd = 1:length(direction)
    
    % Select the proper indices based on cast direction (up or down)
    if strcmp(direction{dd},'up')
        theind = ind_up;
    elseif strcmp(direction{dd},'down')
        theind = ind_down;
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Process Oxygen data
    %%%%%%%%%%%%%%%%%%%%%
    
    % Step 0 - extract variables from RSK file
    % Pressure
    OXY_pres = RSK.data.values(theind,3);
    % Time
    % Correct for date issue on Section 16
    if Sec_num == 16
        OXY_time = RSK.data.tstamp(theind)+6724.46374999;
    else
        OXY_time = RSK.data.tstamp(theind);
    end
    % Oxygen Saturation
    OXY_O2sat = RSK.data.values(theind,4);
    
    %% STEP 1 - Glitch removal. aka "Zero-order hold correction"
    % Must test different channels separately.
    %     dpdt = diff(OXY_pres); rmv = find(dpdt == 0)+1; OXY_pres(rmv) = NaN;
    %     dodt = diff(OXY_O2sat); rmv = find(dodt == 0)+1; OXY_O2sat(rmv) = NaN;
    %     dtdt = diff(OXY_time); rmv = find(dtdt == 0)+1; OXY_time(rmv) = NaN;
    dpdt = diff(OXY_pres);dodt = diff(OXY_O2sat);dtdt = diff(OXY_time);
    ind = find(dpdt == 0 | dodt == 0 | dtdt == 0);
    OXY_pres(ind+1) = NaN;
    OXY_O2sat(ind+1) = NaN;
    clear dpdt dodt dtdt rmv ind
    
    %% STEP 2 - Apply lag correction
    
    % First, compute the lag due to physical separation of the sensors
    % Fall rate
    ws = cat(1,0,diff(OXY_pres)./(diff(OXY_time)*86400));
    % Lag in number of scans (sensor separation = 0.48 m, 8Hz sampling)
    lagadv = 0.48./ws*8;
    
    % Second, apply the lag determined by plotting T-O2 plots
    OXY_O2sat(:,1) = interp1(1:length(OXY_O2sat),OXY_O2sat,(1:length(OXY_O2sat))+lagadv'+lag);
    
    %% STEP 3 - Correction for atmospheric pressure
    OXY_SeaPress = OXY_pres-10.1325;
    
    %% STEP 4 - Secondary variables
    % Longitude & Latitude
    % Get boat position at each timestamp
    [lon,lat] = findcoordinates(OXY_time);
    % Only uses first value (i.e. dropping location) for (lon,lat)
    if strcmp(direction{dd},'down')
        theOXYlondown = lon(1);
        theOXYlatdown = lat(1);
        OXY_lon = lon(1)*ones(length(theind),1); clear lon
        OXY_lat = lat(1)*ones(length(theind),1); clear lat
        % Linear interpolation in time between drop point and boat location
        % when put on board.
    elseif strcmp(direction{dd},'up')
        OXY_lon = interp1([OXY_time(1) OXY_time(end)],[theOXYlondown lon(end)],OXY_time);
        OXY_lat = interp1([OXY_time(1) OXY_time(end)],[theOXYlatdown lat(end)],OXY_time);
        clear theOXYlondown theOXYlatdown
    end
    
    % Oxygen concentration in mL/L and umol/L
    % Read corrected T & S from netcdf file saved for CTD data
    ctdfilename = [filepath 'EcoCTD_ctd/EcoCTD_ctd_S',num2str(Sec_num,'%02d'),...
        '_P',num2str(P_num,'%03d'),'_',direction{dd},'.nc'];
    
    OXY_O2(:,1) = OXY_O2sat/100.*gsw_O2sol(ncread(ctdfilename,'SA'),ncread(ctdfilename,'CT'),ncread(ctdfilename,'seaPress'),ncread(ctdfilename,'lon'),ncread(ctdfilename,'lat'));
    
    %% Step 5 - Apply cross-calibration factors
    OXY_O2sat = SF*OXY_O2sat+offset;
    
    % Recompute secondary variables
    OXY_O2(:,1) = OXY_O2sat/100.*gsw_O2sol(ncread(ctdfilename,'SA'),ncread(ctdfilename,'CT'),ncread(ctdfilename,'seaPress'),ncread(ctdfilename,'lon'),ncread(ctdfilename,'lat'));
    clear ctdfilename
    
    %% STEP 6 - Applies a 16 point median filter to oxygen (~2s)
    %OXY_O2sat = medfilt1(OXY_O2sat,16);
    %OXY_O2mL = medfilt1(OXY_O2mL,16);
    %OXY_O2umol = medfilt1(OXY_O2umol,16);
    
    %% Step 7 - Set the value flag
    OXY_O2sat_flag = ones(length(OXY_time),1);
    OXY_O2_flag = ones(length(OXY_time),1);
    
    % Put a zero Flag is drop rate is too slow
    OXY_O2sat_flag(abs(ws)<1) = 0;
    OXY_O2_flag(abs(ws)<1) = 0;
    % Zero flag if less than 10 m on the way down, 20 m on the way up
    if strcmp(direction,'up')
        OXY_O2sat_flag = -1;
        OXY_O2_flag = -1;
    elseif strcmp(direction,'down')
        OXY_O2sat_flag(OXY_SeaPress<10) = 0;
        OXY_O2_flag(OXY_SeaPress<10) = 0;
    end
    
    % Cap was left on for those casts
    if P_num>=37 && P_num<=45
        OXY_O2sat = NaN*ones(length(OXY_O2sat),1);
        OXY_O2 = NaN*ones(length(OXY_O2sat),1);
        OXY_O2sat_flag = -1*ones(length(OXY_O2sat),1);
        OXY_O2_flag = -1*ones(length(OXY_O2sat),1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Generate NetCDF file
    %%%%%%%%%%%%%%%%%%%%%
    
    % Generate filename
    filename = [filepath '/EcoCTD_oxy/EcoCTD_oxy_S',num2str(Sec_num,'%02d'),...
        '_P',num2str(P_num,'%03d'),'_',direction{dd},'.nc'];
    
    % Create the NetCDF file
    ncid = netcdf.create(filename,'NOCLOBBER');
    
    % Define dimensions
    time_dimID = netcdf.defDim(ncid,'Time',length(theind));
    
    % Define Global Attributes
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'title','Level-1 processed data from Oxygen instrument on EcoCTD package');
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
    netcdf.putAtt(ncid,Pid,'actual_range',[min(OXY_pres(:)) max(OXY_pres(:))]);
    netcdf.putAtt(ncid,Pid,'missing_value',-999);
    netcdf.putAtt(ncid,Pid,'origin','Measured');
    
    tid = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,tid,'long_name','Time');
    netcdf.putAtt(ncid,tid,'units','days');
    netcdf.putAtt(ncid,tid,'axis','T');
    netcdf.putAtt(ncid,tid,'valid_range',[min(OXY_time(:)) max(OXY_time(:))]);
    netcdf.putAtt(ncid,tid,'actual_range',[min(OXY_time(:)) max(OXY_time(:))]);
    netcdf.putAtt(ncid,tid,'missing_value',-999);
    netcdf.putAtt(ncid,tid,'origin','Measured');
    netcdf.putAtt(ncid,tid,'notes','days since [0000 0 0 0 0 0]');
    
    lonid = netcdf.defVar(ncid,'lon','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,lonid,'long_name','Longitude');
    netcdf.putAtt(ncid,lonid,'standard_name','longitude');
    netcdf.putAtt(ncid,lonid,'units','degree_east');
    netcdf.putAtt(ncid,lonid,'valid_range',[-180 180]);
    netcdf.putAtt(ncid,lonid,'actual_range',[min(OXY_lon(:)) max(OXY_lon(:))]);
    netcdf.putAtt(ncid,lonid,'missing_value',-999);
    netcdf.putAtt(ncid,lonid,'origin','Computed');
    
    latid = netcdf.defVar(ncid,'lat','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,latid,'long_name','Latitude');
    netcdf.putAtt(ncid,latid,'standard_name','latitude');
    netcdf.putAtt(ncid,latid,'units','degree_north');
    netcdf.putAtt(ncid,latid,'valid_range',[-90 90]);
    netcdf.putAtt(ncid,latid,'actual_range',[min(OXY_lat(:)) max(OXY_lat(:))]);
    netcdf.putAtt(ncid,latid,'missing_value',-999);
    netcdf.putAtt(ncid,latid,'origin','Computed');
    
    sPid = netcdf.defVar(ncid,'seaPress','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,sPid,'long_name','Sea Pressure');
    netcdf.putAtt(ncid,sPid,'standard_name','sea_water_pressure_due_to_sea_water');
    netcdf.putAtt(ncid,sPid,'units','dbar');
    netcdf.putAtt(ncid,sPid,'valid_range',[0 600]);
    netcdf.putAtt(ncid,sPid,'actual_range',[min(OXY_SeaPress(:)) max(OXY_SeaPress(:))]);
    netcdf.putAtt(ncid,sPid,'missing_value',-999);
    netcdf.putAtt(ncid,sPid,'origin','Computed');
    netcdf.putAtt(ncid,sPid,'notes','Computed as P-10.1325 dbar');
    
    O2satid = netcdf.defVar(ncid,'O2_sat','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,O2satid,'long_name','Oxygen Saturation');
    netcdf.putAtt(ncid,O2satid,'standard_name','fractional_saturation_of_oxygen_in_sea_water');
    netcdf.putAtt(ncid,O2satid,'units','percent');
    netcdf.putAtt(ncid,O2satid,'valid_range',[0 100]);
    netcdf.putAtt(ncid,O2satid,'actual_range',[min(OXY_O2sat(:)) max(OXY_O2sat(:))]);
    netcdf.putAtt(ncid,O2satid,'missing_value',-999);
    netcdf.putAtt(ncid,O2satid,'ancillary_variables','O2_sat_qc');
    netcdf.putAtt(ncid,O2satid,'origin','Computed');
    netcdf.putAtt(ncid,O2satid,'notes',['Calibrated against Ship CTD, scale factor = ',num2str(SF),', offset = ',num2str(offset)]);
    
    O2satflagid = netcdf.defVar(ncid,'O2_sat_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,O2satflagid,'long_name','Oxygen Saturation Quality');
    netcdf.putAtt(ncid,O2satflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,O2satflagid,'units','1');
    netcdf.putAtt(ncid,O2satflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,O2satflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,O2satflagid,'missing_value',-999);
    netcdf.putAtt(ncid,O2satflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,O2satflagid,'flag_meanings','bad questionable good');
    
    O2umolid = netcdf.defVar(ncid,'O2_umolkg','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,O2umolid,'long_name','Oxygen concentration in micromolar per kilogram');
    netcdf.putAtt(ncid,O2umolid,'standard_name','moles_of_oxygen_per_unit_mass_in_sea_water');
    netcdf.putAtt(ncid,O2umolid,'units','micromole kilogram-1');
    netcdf.putAtt(ncid,O2umolid,'valid_range',[0 max(OXY_O2(:))]);
    netcdf.putAtt(ncid,O2umolid,'actual_range',[min(OXY_O2(:)) max(OXY_O2(:))]);
    netcdf.putAtt(ncid,O2umolid,'missing_value',-999);
    netcdf.putAtt(ncid,O2umolid,'ancillary_variables','O2_umolkg_qc');
    netcdf.putAtt(ncid,O2umolid,'origin','Computed');
    netcdf.putAtt(ncid,O2umolid,'notes','Computed using Gibbs SeaWater routine gsw_O2sol');
    
    O2umolflagid = netcdf.defVar(ncid,'O2_umolkg_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,O2umolflagid,'long_name','Oxygen concentration in micromolar per kilogram Quality');
    netcdf.putAtt(ncid,O2umolflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,O2umolflagid,'units','1');
    netcdf.putAtt(ncid,O2umolflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,O2umolflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,O2umolflagid,'missing_value',-999);
    netcdf.putAtt(ncid,O2umolflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,O2umolflagid,'flag_meanings','bad questionable good');
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Write Data to Netcdf File
    %%%%%%%%%%%%%%%%%%%%%
    
    % Absolute Pressure
    netcdf.putVar(ncid,Pid,OXY_pres);
    
    % Time
    netcdf.putVar(ncid,tid,OXY_time);
    
    % Sea Pressure
    netcdf.putVar(ncid,sPid,OXY_SeaPress);
    
    % Longitude & Latitude
    netcdf.putVar(ncid,lonid,OXY_lon);
    netcdf.putVar(ncid,latid,OXY_lat);
    
    % Oxygen Saturation
    netcdf.putVar(ncid,O2satid,OXY_O2sat);
    
    % Oxygen Saturation
    netcdf.putVar(ncid,O2satflagid,OXY_O2sat_flag);
    
    % Oxygen concentration
    netcdf.putVar(ncid,O2umolid,OXY_O2);
    
    % Oxygen concentration flag
    netcdf.putVar(ncid,O2umolflagid,OXY_O2_flag);
    
    netcdf.close(ncid);
    clear filename theind ncid *ID *id OXY*
    
end
