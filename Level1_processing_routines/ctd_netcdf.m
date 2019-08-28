%% This code generates the Level-1 NetCDF files for the EcoCTD.
% It creates a file for downcasts and upcasts for the variables
% measured -- and computed from -- the CTD probe.
%
%
filepath = './'
%% Processing steps are
%       1- Glitch removal.
%           This was attributed to "Analog to Digital (A2D)
%           calibration, which repeats previously logged values. Exact
%           repetition of a datapoint in pressure, conductivity, or temperature
%           are removed from the data.
%
%       2- Apply lag correction
%           a specific lag value was determined by Mara Freilich to minimize
%           the spiking in the data due to lag between temperature and salinity
%           probes at sharp gradients. Constant necessary for processing
%           is:
lag = 0.5;
%
%       3- Apply offset from Guard.
%           CTD was calibrated before equipped with a plastic guard that
%           introduced an offset in the conductivity. The offset is
%           computed based on CTD calibration profiles with the EcoCTD
%           mounted on the rosette, and is estimated to:
temp_SF = 0.98125;
temp_offset = 0.32092;
cond_offset = 0.04361;
%
%       4- Correction for atmospheric pressure (minus 10.1325dbar)
%
%       5- Secondary variables are computed using the TEOS-10 Gibbs
%           Seawater package. These variables are:
%               - Sea pressure (SP = P - 10.1325)
%               - Conservative Temperature (gsw_CT_from_t)
%               - Practical Salinity (gsw_SP_from_C)
%               - Absolute Salinity (gsw_SA_from_SP)
%               - In-situ density (gsw_rho)
%               - Longitude and Latitude. On downcast: boat position at
%               drop point. On upcast: linear interpolation in time between
%               drop location at beginning of upcast, and boat location at
%               end of upcast.
%
%       6 - Applies a 3 point median filter to smooth salinity

direction = {'down','up'};

for dd = 1:length(direction)
    
    % Select the proper indices based on cast direction (up or down)
    if strcmp(direction{dd},'up')
        theind = ind_up;
    elseif strcmp(direction{dd},'down')
        theind = ind_down;
    end
        
    %%%%%%%%%%%%%%%%%%%%%
    %% Process CTD data
    %%%%%%%%%%%%%%%%%%%%%
    
    %% STEP 0 - extract variables from RSK file
    % Pressure
    CTD_pres = RSK.data.values(theind,3);
    % Time
    % Correct for date issue on Section 16
    if Sec_num == 16
        CTD_time = RSK.data.tstamp(theind)+6724.46374999;
    else
        CTD_time = RSK.data.tstamp(theind);
    end
    % Temperature
    CTD_temp = RSK.data.values(theind,2);
    % Conductivity
    CTD_cond = RSK.data.values(theind,1);
    
    %% STEP 1 - Glitch removal. aka "Zero-order hold correction"
    % Must test different channels separately.
    %dpdt = diff(CTD_pres); rmv = find(dpdt == 0)+1; CTD_pres(rmv) = NaN;
    %dcdt = diff(CTD_cond); rmv = find(dcdt == 0)+1; CTD_cond(rmv) = NaN;
    %dtdt = diff(CTD_temp); rmv = find(dtdt == 0)+1; CTD_temp(rmv) = NaN;
    % Change of strategy - filter applied if any of the variables repeat
    dpdt = diff(CTD_pres);dcdt = diff(CTD_cond); dtdt = diff(CTD_temp);
    ind = find(dpdt == 0 | dcdt == 0 | dtdt == 0);
    CTD_pres(ind+1) = NaN;
    CTD_cond(ind+1) = NaN;
    CTD_temp(ind+1) = NaN; 
    clear dpdt dcdt dtdt rmv ind
    
    %% STEP 2 - Apply lag correction. aka "Depsiking"
    CTD_cond(:,1) = interp1(1:length(CTD_cond),CTD_cond,(1:length(CTD_cond))-lag);
    
    %% STEP 3 - Apply cross-correlation
    CTD_temp = temp_SF*CTD_temp+temp_offset;
    CTD_cond = CTD_cond+cond_offset;
    
    %% STEP 4 - Correction for atmospheric pressure
    CTD_SeaPress = CTD_pres-10.1325;
    
    %% STEP 5 - Computed variables
    % Longitude & Latitude
    % Get boat position at each timestamp
    [lon,lat] = findcoordinates(CTD_time);
    % Only uses first value (i.e. dropping location) for (lon,lat)
    if strcmp(direction{dd},'down')
        thectdlondown = lon(1);
        thectdlatdown = lat(1);
        if isnan(thectdlondown)==1 || isnan(thectdlatdown)==1
            error('No longitude/latitude was found for this profile!')
        end
        CTD_lon = lon(1)*ones(length(theind),1); clear lon
        CTD_lat = lat(1)*ones(length(theind),1); clear lat
        % Linear interpolation in time between drop point and boat location
        % when put on board.
    elseif strcmp(direction{dd},'up')
        CTD_lon = interp1([CTD_time(1) CTD_time(end)],[thectdlondown lon(end)],CTD_time);
        CTD_lat = interp1([CTD_time(1) CTD_time(end)],[thectdlatdown lat(end)],CTD_time);
        clear thectdlondown thectdlatdown
    end
    
    % Practical Salinity
    CTD_SP(:,1) = gsw_SP_from_C(CTD_cond,CTD_temp,CTD_SeaPress);
    
    % Absolute Salinity
    [CTD_SA(:,1), ~] = gsw_SA_from_SP(CTD_SP,CTD_SeaPress,CTD_lon,CTD_lat);
    
    % Conservative Temperature
    CTD_CT = gsw_CT_from_t(CTD_SA,CTD_temp,CTD_SeaPress);
    
    %% STEP 6 - Applies a 3 point median filter to salinity
    CTD_SP = medfilt1(CTD_SP,3);
    CTD_SA = medfilt1(CTD_SA,3);
    
    % In-situ density
    CTD_rho = gsw_rho(CTD_SA,CTD_CT,CTD_SeaPress);
    
    %% STEP 7 - Applies flags
    
    CTD_temp_flag = ones(length(CTD_time),1);
    CTD_cond_flag = ones(length(CTD_time),1);
    CTD_CT_flag = ones(length(CTD_time),1);
    CTD_SP_flag = ones(length(CTD_time),1);
    CTD_SA_flag = ones(length(CTD_time),1);
    CTD_rho_flag = ones(length(CTD_time),1);
    
    % Zero flag if less than 10 m on the way down, whole profile on way up
    if strcmp(direction{dd},'up')
        CTD_temp_flag(:) = -1;
        CTD_cond_flag(:) = -1;
        CTD_CT_flag(:) = -1;
        CTD_SP_flag(:) = -1;
        CTD_SA_flag(:) = -1;
        CTD_rho_flag(:) = -1;
    elseif strcmp(direction{dd},'down')
        CTD_temp_flag(CTD_SeaPress<10) = 0;
        CTD_cond_flag(CTD_SeaPress<10) = 0;
        CTD_CT_flag(CTD_SeaPress<10) = 0;
        CTD_SP_flag(CTD_SeaPress<10) = 0;
        CTD_SA_flag(CTD_SeaPress<10) = 0;
        CTD_rho_flag(CTD_SeaPress<10) = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%
    %% Create Netcdf File
    %%%%%%%%%%%%%%%%%%%%%
    
        % Generate filename
    filename = [filepath 'EcoCTD_ctd/EcoCTD_ctd_S',num2str(Sec_num,'%02d'),...
        '_P',num2str(P_num,'%03d'),'_',direction{dd},'.nc'];

    % Create the NetCDF file
    ncid = netcdf.create(filename,'NOCLOBBER');
    
    % Define dimensions
    time_dimID = netcdf.defDim(ncid,'Time',length(theind));
    
    % Define Global Attributes
    varid = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid,varid,'title','Level-1 processed data from CTD instrument on EcoCTD package');
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
    netcdf.putAtt(ncid,Pid,'actual_range',[min(CTD_pres(:)) max(CTD_pres(:))]);
    netcdf.putAtt(ncid,Pid,'missing_value',-999);
    netcdf.putAtt(ncid,Pid,'origin','Measured');
    
    tid = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,tid,'long_name','Time');
    netcdf.putAtt(ncid,tid,'units','days');
    netcdf.putAtt(ncid,tid,'axis','T');
    netcdf.putAtt(ncid,tid,'valid_range',[min(CTD_time(:)) max(CTD_time(:))]);
    netcdf.putAtt(ncid,tid,'actual_range',[min(CTD_time(:)) max(CTD_time(:))]);
    netcdf.putAtt(ncid,tid,'missing_value',-999);
    netcdf.putAtt(ncid,tid,'origin','Measured');
    netcdf.putAtt(ncid,tid,'notes','days since [0000 0 0 0 0 0]');
    
    lonid = netcdf.defVar(ncid,'lon','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,lonid,'long_name','Longitude');
    netcdf.putAtt(ncid,lonid,'standard_name','longitude');
    netcdf.putAtt(ncid,lonid,'units','degree_east');
    netcdf.putAtt(ncid,lonid,'valid_range',[-180 180]);
    netcdf.putAtt(ncid,lonid,'actual_range',[min(CTD_lon(:)) max(CTD_lon(:))]);
    netcdf.putAtt(ncid,lonid,'missing_value',-999);
    netcdf.putAtt(ncid,lonid,'origin','Computed');
    
    latid = netcdf.defVar(ncid,'lat','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,latid,'long_name','Latitude');
    netcdf.putAtt(ncid,latid,'standard_name','latitude');
    netcdf.putAtt(ncid,latid,'units','degree_north');
    netcdf.putAtt(ncid,latid,'valid_range',[-90 90]);
    netcdf.putAtt(ncid,latid,'actual_range',[min(CTD_lat(:)) max(CTD_lat(:))]);
    netcdf.putAtt(ncid,latid,'missing_value',-999);
    netcdf.putAtt(ncid,latid,'origin','Computed');
    
    Tid = netcdf.defVar(ncid,'T','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,Tid,'long_name','In-situ Temperature');
    netcdf.putAtt(ncid,Tid,'standard_name','sea_water_temperature');
    netcdf.putAtt(ncid,Tid,'units','Celsius');
    netcdf.putAtt(ncid,Tid,'valid_range',[0 40]);
    netcdf.putAtt(ncid,Tid,'actual_range',[min(CTD_temp(:)) max(CTD_temp(:))]);
    netcdf.putAtt(ncid,Tid,'missing_value',-999);
    netcdf.putAtt(ncid,Tid,'origin','Measured');
    
    Tflagid = netcdf.defVar(ncid,'T_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,Tflagid,'long_name','In-situ Temperature Quality');
    netcdf.putAtt(ncid,Tflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,Tflagid,'units','1');
    netcdf.putAtt(ncid,Tflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,Tflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,Tflagid,'missing_value',-999);
    netcdf.putAtt(ncid,Tflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,Tflagid,'flag_meanings','bad questionable good');

    Cid = netcdf.defVar(ncid,'C','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,Cid,'long_name','Conductivity');
    netcdf.putAtt(ncid,Cid,'standard_name','sea_water_electrical_conductivity');
    netcdf.putAtt(ncid,Cid,'units','mS cm-1');
    netcdf.putAtt(ncid,Cid,'valid_range',[min(CTD_cond(:)) max(CTD_cond(:))]);
    netcdf.putAtt(ncid,Cid,'actual_range',[min(CTD_cond(:)) max(CTD_cond(:))]);
    netcdf.putAtt(ncid,Cid,'missing_value',-999);
    netcdf.putAtt(ncid,Cid,'origin','Measured');

    Cflagid = netcdf.defVar(ncid,'C_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,Cflagid,'long_name','Conductivity Quality');
    netcdf.putAtt(ncid,Cflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,Cflagid,'units','1');
    netcdf.putAtt(ncid,Cflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,Cflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,Cflagid,'missing_value',-999);
    netcdf.putAtt(ncid,Cflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,Cflagid,'flag_meanings','bad questionable good');
    
    sPid = netcdf.defVar(ncid,'seaPress','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,sPid,'long_name','Sea Pressure');
    netcdf.putAtt(ncid,sPid,'standard_name','sea_water_pressure_due_to_sea_water');
    netcdf.putAtt(ncid,sPid,'units','dbar');
    netcdf.putAtt(ncid,sPid,'valid_range',[0 600]);
    netcdf.putAtt(ncid,sPid,'actual_range',[min(CTD_SeaPress(:)) max(CTD_SeaPress(:))]);
    netcdf.putAtt(ncid,sPid,'missing_value',-999);
    netcdf.putAtt(ncid,sPid,'origin','Computed');
    netcdf.putAtt(ncid,sPid,'notes','Computed as P-10.1325 dbar');
    
    CTid = netcdf.defVar(ncid,'CT','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,CTid,'long_name','Conservative Temperature');
    netcdf.putAtt(ncid,CTid,'standard_name','sea_water_conservative_temperature');
    netcdf.putAtt(ncid,CTid,'units','Celsius');
    netcdf.putAtt(ncid,CTid,'valid_range',[0 40]);
    netcdf.putAtt(ncid,CTid,'actual_range',[min(CTD_CT(:)) max(CTD_CT(:))]);
    netcdf.putAtt(ncid,CTid,'missing_value',-999);
    netcdf.putAtt(ncid,CTid,'origin','Computed');
    netcdf.putAtt(ncid,CTid,'notes','Computed using the GibbsSeawater toolbox (gsw_CT_from_t)');

    CTflagid = netcdf.defVar(ncid,'CT_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,CTflagid,'long_name','Conservative Temperature Quality');
    netcdf.putAtt(ncid,CTflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,CTflagid,'units','1');
    netcdf.putAtt(ncid,CTflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,CTflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,CTflagid,'missing_value',-999);
    netcdf.putAtt(ncid,CTflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,CTflagid,'flag_meanings','bad questionable good');
    
    SPid = netcdf.defVar(ncid,'SP','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,SPid,'long_name','Practical Salinity');
    netcdf.putAtt(ncid,SPid,'standard_name','sea_water_practical_salinity');
    netcdf.putAtt(ncid,SPid,'units','1');
    netcdf.putAtt(ncid,SPid,'valid_range',[0 45]);
    netcdf.putAtt(ncid,SPid,'actual_range',[min(CTD_SP(:)) max(CTD_SP(:))]);
    netcdf.putAtt(ncid,SPid,'missing_value',-999);
    netcdf.putAtt(ncid,SPid,'origin','Computed');
    netcdf.putAtt(ncid,SPid,'notes','Computed using the GibbsSeawater toolbox (gsw_SP_from_C)');
    
    SPflagid = netcdf.defVar(ncid,'SP_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,SPflagid,'long_name','Practical Salinity Quality');
    netcdf.putAtt(ncid,SPflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,SPflagid,'units','1');
    netcdf.putAtt(ncid,SPflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,SPflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,SPflagid,'missing_value',-999);
    netcdf.putAtt(ncid,SPflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,SPflagid,'flag_meanings','bad questionable good');
        
    SAid = netcdf.defVar(ncid,'SA','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,SAid,'long_name','Absolute Salinity');
    netcdf.putAtt(ncid,SAid,'standard_name','sea_water_absolute_salinity');
    netcdf.putAtt(ncid,SAid,'units','gram kilogram-1');
    netcdf.putAtt(ncid,SAid,'valid_range',[0 45]);
    netcdf.putAtt(ncid,SAid,'actual_range',[min(CTD_SA(:)) max(CTD_SA(:))]);
    netcdf.putAtt(ncid,SAid,'missing_value',-999);
    netcdf.putAtt(ncid,SAid,'origin','Computed');
    netcdf.putAtt(ncid,SAid,'notes','Computed using the GibbsSeawater toolbox (gsw_SA_from_SP)');

    SAflagid = netcdf.defVar(ncid,'SA_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,SAflagid,'long_name','Absolute Salinity Quality');
    netcdf.putAtt(ncid,SAflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,SAflagid,'units','1');
    netcdf.putAtt(ncid,SAflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,SAflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,SAflagid,'missing_value',-999);
    netcdf.putAtt(ncid,SAflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,SAflagid,'flag_meanings','bad questionable good');
        
    densid = netcdf.defVar(ncid,'rho','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,densid,'long_name','In-situ Density');
    netcdf.putAtt(ncid,densid,'standard_name','sea_water_density');
    netcdf.putAtt(ncid,densid,'units','kilogram meter-3');
    netcdf.putAtt(ncid,densid,'valid_range',[1000 1050]);
    netcdf.putAtt(ncid,densid,'actual_range',[min(CTD_rho(:)) max(CTD_rho(:))]);
    netcdf.putAtt(ncid,densid,'missing_value',-999);
    netcdf.putAtt(ncid,densid,'origin','Computed');
    netcdf.putAtt(ncid,densid,'notes','Computed using the GibbsSeawater toolbox (gsw_rho)');

    densflagid = netcdf.defVar(ncid,'rho_qc','NC_DOUBLE',time_dimID);
    netcdf.putAtt(ncid,densflagid,'long_name','In-situ Density Quality');
    netcdf.putAtt(ncid,densflagid,'standard_name','status_flag');
    netcdf.putAtt(ncid,densflagid,'units','1');
    netcdf.putAtt(ncid,densflagid,'valid_range',[-1 1]);
    netcdf.putAtt(ncid,densflagid,'actual_range',[-1 1]);
    netcdf.putAtt(ncid,densflagid,'missing_value',-999);
    netcdf.putAtt(ncid,densflagid,'flag_values',[-1 0 1]);
    netcdf.putAtt(ncid,densflagid,'flag_meanings','bad questionable good');        
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid);
    
    %%%%%%%%%%%%%%%%%%%%%
    %% Write Data to Netcdf File
    %%%%%%%%%%%%%%%%%%%%%
    
    % Absolute Pressure
    netcdf.putVar(ncid,Pid,CTD_pres);
    
    % Time
    netcdf.putVar(ncid,tid,CTD_time);
    
    % Longitude & Latitude
    netcdf.putVar(ncid,lonid,CTD_lon);
    netcdf.putVar(ncid,latid,CTD_lat);
    
    % in-situ Temperature
    netcdf.putVar(ncid,Tid,CTD_temp);
    
    % in-situ Temperature flag
    netcdf.putVar(ncid,Tflagid,CTD_temp_flag);    
    
    % Conductivity
    netcdf.putVar(ncid,Cid,CTD_cond);
    
    % Conductivity flag
    netcdf.putVar(ncid,Cflagid,CTD_cond_flag);    
    
    % Sea Pressure
    netcdf.putVar(ncid,sPid,CTD_SeaPress);
    
    % Conservative temperature
    netcdf.putVar(ncid,CTid,CTD_CT);
    
    % Conservative temperature flag
    netcdf.putVar(ncid,CTflagid,CTD_CT_flag);
    
    % Practical Salinity
    netcdf.putVar(ncid,SPid,CTD_SP);
    
    % Practical Salinity flag
    netcdf.putVar(ncid,SPflagid,CTD_SP_flag);    
    
    % Absolute Salinity
    netcdf.putVar(ncid,SAid,CTD_SA);
    
    % Absolute Salinity flag
    netcdf.putVar(ncid,SAflagid,CTD_SA_flag);    
    
    % In-situ density
    netcdf.putVar(ncid,densid,CTD_rho);
    
    % In-situ density flag
    netcdf.putVar(ncid,densflagid,CTD_rho_flag);    
    
    netcdf.close(ncid);
    clear filename theind ncid *ID *id CTD*
    
end; clear direction dd
