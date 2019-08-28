%% This code generates the Level-3 NetCDF files for the merged EcoCTD and UCTD.
clear
%%%%%%%%%%%%%%%%%%%%%
%% Open Level 1 Data from Netcdf File
%%%%%%%%%%%%%%%%%%%%%

% Define depth grid
P = 0:.5:300;
% Path to ecoCTD files
path = './'

%% EcoCTD data

%%%%%%%%%%%%%%%%%%%%%
% EcoCTD CTD data
%%%%%%%%%%%%%%%%%%%%%
files = dir([path 'EcoCTD_ctd/*down*']);

for ii = 1:length(files)
    disp(num2str(ii));
    filename = fullfile(files(ii).folder,files(ii).name);
    
    % Get the variables
    eco.P = ncread(filename,'seaPress');
    eco.time = ncread(filename,'time');
    eco.lon = ncread(filename,'lon');
    eco.lat = ncread(filename,'lat');
    eco.T = ncread(filename,'T');
    eco.T_flag = ncread(filename,'T_qc');
    eco.CT = ncread(filename,'CT');
    eco.CT_flag = ncread(filename,'CT_qc');
    eco.SP = ncread(filename,'SP');
    eco.SP_flag = ncread(filename,'SP_qc');
    eco.SA = ncread(filename,'SA');
    eco.SA_flag = ncread(filename,'SA_qc');
    eco.rho = ncread(filename,'rho');
    eco.rho_flag = ncread(filename,'rho_qc');
    
    % Grid the variables
    fields = fieldnames(eco);
    [eco.P,IA,~] = unique(eco.P);
    for jj = 2:length(fields)
        eval(['eco.',fields{jj},' = eco.',fields{jj},'(IA);'])
        eval(['good = find(isnan(eco.',fields{1},') == 0 & isnan(eco.',fields{jj},') == 0 );'])
        if isempty(good)
            eval(['ctd.',fields{jj},'(:,ii) = NaN*ones(length(P),1);'])
        else
            eval(['ctd.',fields{jj},'(:,ii) = interp1(eco.',fields{1},'(good),eco.',fields{jj},'(good),P,''spline'');'])
        end
    end; clear jj
    
end; clear ii eco

% Replace NaNs by FillValue
fields = fieldnames(ctd);
for jj = 2:length(fields)
    % -1 flag for flag variables, -999 value for the rest
    if strcmp(fields{jj}(end),'g')
        eval(['ctd.',fields{jj},'(isnan(ctd.',fields{jj},')==1) = -1;'])
    else
        eval(['ctd.',fields{jj},'(isnan(ctd.',fields{jj},')==1) = -999;'])
    end
end; clear jj

%%%%%%%%%%%%%%%%%%%%%
% EcoCTD OXY data
%%%%%%%%%%%%%%%%%%%%%
files = dir([path 'EcoCTD_oxy/*down*']);

for ii = 1:length(files)
    disp(num2str(ii));
    filename = fullfile(files(ii).folder,files(ii).name);
    
    % Get the variables
    eco.P = ncread(filename,'sP');
    eco.time = ncread(filename,'time');
    eco.lon = ncread(filename,'lon');
    eco.lat = ncread(filename,'lat');
    eco.O2_sat = ncread(filename,'O2_sat');
    eco.O2_sat_flag = ncread(filename,'O2_sat_qc');
    eco.O2_umolkg = ncread(filename,'O2_umolkg');
    eco.O2_umolkg_flag = ncread(filename,'O2_umolkg_qc');
    
    % Grid the variables
    fields = fieldnames(eco);
    [eco.P,IA,~] = unique(eco.P);
    for jj = 2:length(fields)
        eval(['eco.',fields{jj},' = eco.',fields{jj},'(IA);'])
        eval(['good = find(isnan(eco.',fields{1},') == 0 & isnan(eco.',fields{jj},') == 0 );'])
        if isempty(good)
            eval(['oxy.',fields{jj},'(:,ii) = NaN*ones(length(P),1);'])
        else
            eval(['oxy.',fields{jj},'(:,ii) = interp1(eco.',fields{1},'(good),eco.',fields{jj},'(good),P,''spline'');'])
        end
    end; clear jj
    
end; clear ii eco

% Replace NaNs by FillValue
fields = fieldnames(oxy);
for jj = 2:length(fields)
    % -1 flag for flag variables, -999 value for the rest
    if strcmp(fields{jj}(end-1:end),'ag')
        eval(['oxy.',fields{jj},'(isnan(oxy.',fields{jj},')==1) = -1;'])
    else
        eval(['oxy.',fields{jj},'(isnan(oxy.',fields{jj},')==1) = -999;'])
    end
end; clear jj

%%%%%%%%%%%%%%%%%%%%%
% EcoCTD FLS data
%%%%%%%%%%%%%%%%%%%%%
files = dir([path 'EcoCTD_fls/*down*']);

for ii = 1:length(files)
    disp(num2str(ii));
    filename = fullfile(files(ii).folder,files(ii).name);
    
    % Get the variables
    eco.P = ncread(filename,'sP');
    eco.time = ncread(filename,'time');
    eco.lon = ncread(filename,'lon');
    eco.lat = ncread(filename,'lat');
    eco.bb470 = ncread(filename,'bb470');
    eco.bb470_flag = ncread(filename,'bb470_qc');
    eco.bb700 = ncread(filename,'bb700');
    eco.bb700_flag = ncread(filename,'bb700_qc');
    eco.chl_cal = ncread(filename,'chl_cal');
    eco.chl_cal_flag = ncread(filename,'chl_cal_qc');
    
    % Despike backscatter variables
    % Running minimum filter
    eco.bb470 = movmin(eco.bb470,7);
    eco.bb700 = movmin(eco.bb700,7);
    % Running maximum filter
    eco.bb470 = movmax(eco.bb470,7);
    eco.bb700 = movmax(eco.bb700,7);
    
    % Grid the variables
    fields = fieldnames(eco);
    [eco.P,IA,~] = unique(eco.P);
    for jj = 2:length(fields)
        eval(['eco.',fields{jj},' = eco.',fields{jj},'(IA);'])
        eval(['good = find(isnan(eco.',fields{1},') == 0 & isnan(eco.',fields{jj},') == 0 );'])
        if isempty(good)
            eval(['fls.',fields{jj},'(:,ii) = NaN*ones(length(P),1);'])
        else
            eval(['fls.',fields{jj},'(:,ii) = interp1(eco.',fields{1},'(good),eco.',fields{jj},'(good),P,''spline'');'])
        end
    end; clear jj
    
end; clear ii eco

% Replace NaNs by FillValue
fields = fieldnames(fls);
for jj = 2:length(fields)
    % -1 flag for flag variables, -999 value for the rest
    if strcmp(fields{jj}(end),'g')
        eval(['fls.',fields{jj},'(isnan(fls.',fields{jj},')==1) = -1;'])
    else
        eval(['fls.',fields{jj},'(isnan(fls.',fields{jj},')==1) = -999;'])
    end
end; clear jj

% Replace NaNs by FillValue
fields = fieldnames(ctd2);
for jj = 2:length(fields)
    % -1 flag for flag variables, -999 value for the rest
    if strcmp(fields{jj}(end),'g')
        eval(['ctd2.',fields{jj},'(isnan(ctd2.',fields{jj},')==1) = -1;'])
    else
        eval(['ctd2.',fields{jj},'(isnan(ctd2.',fields{jj},')==1) = -999;'])
    end
end; clear jj

%%%%%%%%%%%%%%%%%%%%%
%% Organize Data
%%%%%%%%%%%%%%%%%%%%%

% sort chronologically
blop = cat(2,ctd.time,ctd2.time);
ID_sensor = cat(2,ones(size(ctd.time)),2*ones(size(ctd2.time)));
ind = find(isnan(sum(blop,2))==0,1,'first');
[~,I] = sort(blop(ind,:));

% Sort time
CTD_time = blop(:,I);

% Sort time
ID_sensor = ID_sensor(:,I);

% Sort longitude
blop = cat(2,ctd.lon,ctd2.lon);
CTD_lon = blop(:,I);

% Sort latitude
blop = cat(2,ctd.lat,ctd2.lat);
CTD_lat = blop(:,I);

% Sort Temperature
blop = cat(2,ctd.T,ctd2.T);
CTD_T = blop(:,I);

% Sort Temperature flag
blop = cat(2,ctd.T_flag,ctd2.T_flag);
CTD_T_flag = blop(:,I);

% Sort conservative temperature
blop = cat(2,ctd.CT,ctd2.CT);
CTD_CT = blop(:,I);

% Sort conservative temperature flag
blop = cat(2,ctd.CT_flag,ctd2.CT_flag);
CTD_CT_flag = blop(:,I);

% Sort practical salinity
blop = cat(2,ctd.SP,ctd2.SP);
CTD_SP = blop(:,I);

% Sort practical salinity flag
blop = cat(2,ctd.SP_flag,ctd2.SP_flag);
CTD_SP_flag = blop(:,I);

% Sort absolute salinity
blop = cat(2,ctd.SA,ctd2.SA);
CTD_SA = blop(:,I);

% Sort absolute salinity flag
blop = cat(2,ctd.SA_flag,ctd2.SA_flag);
CTD_SA_flag = blop(:,I);

% Sort density
blop = cat(2,ctd.rho,ctd2.rho);
CTD_rho = blop(:,I);

% Sort density flag
blop = cat(2,ctd.rho_flag,ctd2.rho_flag);
CTD_rho_flag = blop(:,I);

% Sort chlorophyll. Concatenate with missing value as no CHL on UCTD
blop = cat(2,fls.chl_cal,-999*ones(size(ctd2.time)));
FLS_chl = blop(:,I);

% Sort chlorophyll_flag. Concatenate with -1 flag as no CHL on UCTD
blop = cat(2,fls.chl_cal_flag,-1*ones(size(ctd2.time)));
FLS_chl_flag = blop(:,I);

% Sort backscatter at 470nm. Concatenate with empty matrix as no backscatter on UCTD
blop = cat(2,fls.bb470,-999*ones(size(ctd2.time)));
FLS_bb470 = blop(:,I);

% backscatter 470 flag. Same as CHL flag
blop = cat(2,fls.bb470_flag,-1*ones(size(ctd2.time)));
FLS_bb470_flag = blop(:,I);

% Sort backscatter at 700nm. Concatenate with empty matrix as no backscatter on UCTD
blop = cat(2,fls.bb700,-999*ones(size(ctd2.time)));
FLS_bb700 = blop(:,I);

% backscatter 700 flag. Same as CHL flag
blop = cat(2,fls.bb700_flag,-1*ones(size(ctd2.time)));
FLS_bb700_flag = blop(:,I);

% Sort oxygen saturation. Concatenate with empty matrix as no oxygen on UCTD
blop = cat(2,oxy.O2_sat,-999*ones(size(ctd2.time)));
OXY_O2_sat = blop(:,I);

% Sort oxygen saturation flag. Concatenate with empty matrix as no oxygen on UCTD
blop = cat(2,oxy.O2_sat_flag,-1*ones(size(ctd2.time)));
OXY_O2_sat_flag = blop(:,I);

% Sort oxygen concentration. Concatenate with empty matrix as no oxygen on UCTD
blop = cat(2,oxy.O2_umolkg,-999*ones(size(ctd2.time)));
OXY_O2_umolkg = blop(:,I);

% oxygen concentration flag. Same as O2 saturation
blop = cat(2,oxy.O2_umolkg_flag,-1*ones(size(ctd2.time)));
OXY_O2_umolkg_flag = blop(:,I);

%%%%%%%%%%%%%%%%%%%%%
%% Generate NetCDF file
%%%%%%%%%%%%%%%%%%%%%

% Generate filename
filename = [path 'Level3/EcoCTD_UCTD_merged.nc'];

% Create the NetCDF file
ncid = netcdf.create(filename,'NOCLOBBER');

% Define dimensions
pnum_dimID = netcdf.defDim(ncid,'Profile number',size(CTD_time,2));
depth_dimID = netcdf.defDim(ncid,'Depth',size(CTD_time,1));

% Define Global Attributes
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'title','Level-3 product of downcasts from EcoCTD probe.');
netcdf.putAtt(ncid,varid,'institution','Woods Hole Oceanographic Institution');
netcdf.putAtt(ncid,varid,'source','ocean profile observations');
netcdf.putAtt(ncid,varid,'history',[datestr(now),' - File generated by Dr. M. Dever']);
netcdf.putAtt(ncid,varid,'references','CALYPSO Data Report');
netcdf.putAtt(ncid,varid,'external_variables','');
netcdf.putAtt(ncid,varid,'Conventions','CF-1.7');
netcdf.putAtt(ncid,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid,varid,'Comment','');

% Define variables and their attributes
pnumid = netcdf.defVar(ncid,'Profile_number','NC_DOUBLE',pnum_dimID);
netcdf.putAtt(ncid,pnumid,'long_name','Profile number');
netcdf.putAtt(ncid,pnumid,'units','1');
netcdf.putAtt(ncid,pnumid,'valid_range',[1 size(CTD_time,2)]);
netcdf.putAtt(ncid,pnumid,'actual_range',[1 size(CTD_time,2)]);
netcdf.putAtt(ncid,pnumid,'missing_value',-999);

depthid = netcdf.defVar(ncid,'Depth','NC_DOUBLE',depth_dimID);
netcdf.putAtt(ncid,depthid,'long_name','Depth grid');
netcdf.putAtt(ncid,depthid,'units','meters');
netcdf.putAtt(ncid,depthid,'positive','down');
netcdf.putAtt(ncid,depthid,'axis','Z');
netcdf.putAtt(ncid,depthid,'valid_range',[min(P) max(P)]);
netcdf.putAtt(ncid,depthid,'actual_range',[min(P) max(P)]);

sensorid = netcdf.defVar(ncid,'Instrument_ID','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,sensorid,'long_name','ID number of instrument used');
netcdf.putAtt(ncid,sensorid,'standard_name','platform_id');
netcdf.putAtt(ncid,sensorid,'units','1');
netcdf.putAtt(ncid,sensorid,'valid_range',[1 2]);
netcdf.putAtt(ncid,sensorid,'actual_range',[1 2]);
netcdf.putAtt(ncid,sensorid,'missing_value',-999);
netcdf.putAtt(ncid,sensorid,'notes','1- EcoCTD; 2- UCTD');

tid = netcdf.defVar(ncid,'time','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,tid,'long_name','Time');
netcdf.putAtt(ncid,tid,'units','days');
netcdf.putAtt(ncid,tid,'valid_range',[min(CTD_time(:)) max(CTD_time(:))]);
netcdf.putAtt(ncid,tid,'actual_range',[min(CTD_time(:)) max(CTD_time(:))]);
netcdf.putAtt(ncid,tid,'missing_value',-999);
netcdf.putAtt(ncid,tid,'origin','Measured');
netcdf.putAtt(ncid,tid,'notes','days since [0000 0 0 0 0 0]');

lonid = netcdf.defVar(ncid,'lon','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,lonid,'long_name','Longitude');
netcdf.putAtt(ncid,lonid,'standard_name','longitude');
netcdf.putAtt(ncid,lonid,'units','degree_east');
netcdf.putAtt(ncid,lonid,'valid_range',[-180 180]);
netcdf.putAtt(ncid,lonid,'actual_range',[min(CTD_lon(:)) max(CTD_lon(:))]);
netcdf.putAtt(ncid,lonid,'missing_value',-999);
netcdf.putAtt(ncid,lonid,'origin','Computed');

latid = netcdf.defVar(ncid,'lat','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,latid,'long_name','Latitude');
netcdf.putAtt(ncid,latid,'standard_name','latitude');
netcdf.putAtt(ncid,latid,'units','degree_north');
netcdf.putAtt(ncid,latid,'valid_range',[-90 90]);
netcdf.putAtt(ncid,latid,'actual_range',[min(CTD_lat(:)) max(CTD_lat(:))]);
netcdf.putAtt(ncid,latid,'missing_value',-999);
netcdf.putAtt(ncid,latid,'origin','Computed');

Tid = netcdf.defVar(ncid,'T','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,Tid,'long_name','In-situ Temperature');
netcdf.putAtt(ncid,Tid,'standard_name','sea_water_temperature');
netcdf.putAtt(ncid,Tid,'units','Celsius');
netcdf.putAtt(ncid,Tid,'valid_range',[0 40]);
netcdf.putAtt(ncid,Tid,'actual_range',[min(CTD_T(:)) max(CTD_T(:))]);
netcdf.putAtt(ncid,Tid,'missing_value',-999);
netcdf.putAtt(ncid,Tid,'origin','Measured');

Tflagid = netcdf.defVar(ncid,'T_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,Tflagid,'long_name','In-situ Temperature Quality');
netcdf.putAtt(ncid,Tflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,Tflagid,'units','1');
netcdf.putAtt(ncid,Tflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,Tflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,Tflagid,'missing_value',-999);
netcdf.putAtt(ncid,Tflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,Tflagid,'flag_meanings','bad questionable good');

CTid = netcdf.defVar(ncid,'CT','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,CTid,'long_name','Conservative Temperature');
netcdf.putAtt(ncid,CTid,'standard_name','sea_water_conservative_temperature');
netcdf.putAtt(ncid,CTid,'units','Celsius');
netcdf.putAtt(ncid,CTid,'valid_range',[0 40]);
netcdf.putAtt(ncid,CTid,'actual_range',[min(CTD_CT(:)) max(CTD_CT(:))]);
netcdf.putAtt(ncid,CTid,'missing_value',-999);
netcdf.putAtt(ncid,CTid,'origin','Computed');
netcdf.putAtt(ncid,CTid,'notes','Computed using the GibbsSeawater toolbox (gsw_CT_from_t)');

CTflagid = netcdf.defVar(ncid,'CT_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,CTflagid,'long_name','Conservative Temperature Quality');
netcdf.putAtt(ncid,CTflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,CTflagid,'units','1');
netcdf.putAtt(ncid,CTflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,CTflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,CTflagid,'missing_value',-999);
netcdf.putAtt(ncid,CTflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,CTflagid,'flag_meanings','bad questionable good');

SPid = netcdf.defVar(ncid,'SP','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,SPid,'long_name','Practical Salinity');
netcdf.putAtt(ncid,SPid,'standard_name','sea_water_practical_salinity');
netcdf.putAtt(ncid,SPid,'units','1');
netcdf.putAtt(ncid,SPid,'valid_range',[0 45]);
netcdf.putAtt(ncid,SPid,'actual_range',[min(CTD_SP(:)) max(CTD_SP(:))]);
netcdf.putAtt(ncid,SPid,'missing_value',-999);
netcdf.putAtt(ncid,SPid,'origin','Computed');
netcdf.putAtt(ncid,SPid,'notes','Computed using the GibbsSeawater toolbox (gsw_SP_from_C)');

SPflagid = netcdf.defVar(ncid,'SP_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,SPflagid,'long_name','Practical Salinity Quality');
netcdf.putAtt(ncid,SPflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,SPflagid,'units','1');
netcdf.putAtt(ncid,SPflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,SPflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,SPflagid,'missing_value',-999);
netcdf.putAtt(ncid,SPflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,SPflagid,'flag_meanings','bad questionable good');

SAid = netcdf.defVar(ncid,'SA','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,SAid,'long_name','Absolute Salinity');
netcdf.putAtt(ncid,SAid,'standard_name','sea_water_absolute_salinity');
netcdf.putAtt(ncid,SAid,'units','gram kilogram-1');
netcdf.putAtt(ncid,SAid,'valid_range',[0 45]);
netcdf.putAtt(ncid,SAid,'actual_range',[min(CTD_SA(:)) max(CTD_SA(:))]);
netcdf.putAtt(ncid,SAid,'missing_value',-999);
netcdf.putAtt(ncid,SAid,'origin','Computed');
netcdf.putAtt(ncid,SAid,'notes','Computed using the GibbsSeawater toolbox (gsw_SA_from_SP)');

SAflagid = netcdf.defVar(ncid,'SA_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,SAflagid,'long_name','Absolute Salinity Quality');
netcdf.putAtt(ncid,SAflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,SAflagid,'units','1');
netcdf.putAtt(ncid,SAflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,SAflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,SAflagid,'missing_value',-999);
netcdf.putAtt(ncid,SAflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,SAflagid,'flag_meanings','bad questionable good');

densid = netcdf.defVar(ncid,'rho','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,densid,'long_name','In-situ Density');
netcdf.putAtt(ncid,densid,'standard_name','sea_water_density');
netcdf.putAtt(ncid,densid,'units','kilogram meter-3');
netcdf.putAtt(ncid,densid,'valid_range',[1000 1050]);
netcdf.putAtt(ncid,densid,'actual_range',[min(CTD_rho(:)) max(CTD_rho(:))]);
netcdf.putAtt(ncid,densid,'missing_value',-999);
netcdf.putAtt(ncid,densid,'origin','Computed');
netcdf.putAtt(ncid,densid,'notes','Computed using the GibbsSeawater toolbox (gsw_rho)');

densflagid = netcdf.defVar(ncid,'rho_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,densflagid,'long_name','In-situ Density Quality');
netcdf.putAtt(ncid,densflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,densflagid,'units','1');
netcdf.putAtt(ncid,densflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,densflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,densflagid,'missing_value',-999);
netcdf.putAtt(ncid,densflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,densflagid,'flag_meanings','bad questionable good');

bb470id = netcdf.defVar(ncid,'bb470','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,bb470id,'long_name','Despiked Backscatter at 470 nm');
netcdf.putAtt(ncid,bb470id,'units','1');
netcdf.putAtt(ncid,bb470id,'valid_range',[0 max(FLS_bb470(:))]);
netcdf.putAtt(ncid,bb470id,'actual_range',[min(FLS_bb470(:)) max(FLS_bb470(:))]);
netcdf.putAtt(ncid,bb470id,'missing_value',-999);
netcdf.putAtt(ncid,bb470id,'ancillary_variables','bb470_qc');
netcdf.putAtt(ncid,bb470id,'origin','Measured');
netcdf.putAtt(ncid,bb470id,'notes','Non-standard units of [counts]. Profiles have been despiked using minimum and maximum running filter, see Data report');

bb470flagid = netcdf.defVar(ncid,'bb470_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,bb470flagid,'long_name','Despiked Backscatter at 470 nm Quality');
netcdf.putAtt(ncid,bb470flagid,'standard_name','status_flag');
netcdf.putAtt(ncid,bb470flagid,'units','1');
netcdf.putAtt(ncid,bb470flagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,bb470flagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,bb470flagid,'missing_value',-999);
netcdf.putAtt(ncid,bb470flagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,bb470flagid,'flag_meanings','bad questionable good');

bb700id = netcdf.defVar(ncid,'bb700','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,bb700id,'long_name','Despiked Backscatter at 700 nm');
netcdf.putAtt(ncid,bb700id,'units','1');
netcdf.putAtt(ncid,bb700id,'valid_range',[0 max(FLS_bb700(:))]);
netcdf.putAtt(ncid,bb700id,'actual_range',[min(FLS_bb700(:)) max(FLS_bb700(:))]);
netcdf.putAtt(ncid,bb700id,'missing_value',-999);
netcdf.putAtt(ncid,bb700id,'ancillary_variables','bb700_qc');
netcdf.putAtt(ncid,bb700id,'origin','Measured');
netcdf.putAtt(ncid,bb470id,'notes','Non-standard units of [counts]. Profiles have been despiked using minimum and maximum running filter, see Data report');

bb700flagid = netcdf.defVar(ncid,'bb700_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,bb700flagid,'long_name','Despiked Backscatter at 700 nm Quality');
netcdf.putAtt(ncid,bb700flagid,'standard_name','status_flag');
netcdf.putAtt(ncid,bb700flagid,'units','1');
netcdf.putAtt(ncid,bb700flagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,bb700flagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,bb700flagid,'missing_value',-999);
netcdf.putAtt(ncid,bb700flagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,bb700flagid,'flag_meanings','bad questionable good');

chlid = netcdf.defVar(ncid,'chl_cal','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,chlid,'long_name','Calibrated Chlorophyll-Fluorescence');
netcdf.putAtt(ncid,chlid,'units','microgram liter-1');
netcdf.putAtt(ncid,chlid,'valid_range',[0 max(FLS_chl(:))]);
netcdf.putAtt(ncid,chlid,'actual_range',[min(FLS_chl(:)) max(FLS_chl(:))]);
netcdf.putAtt(ncid,chlid,'missing_value',-999);
netcdf.putAtt(ncid,chlid,'ancillary_variables','chl_cal_qc');
netcdf.putAtt(ncid,chlid,'origin','Computed');

chlflagid = netcdf.defVar(ncid,'chl_cal_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,chlflagid,'long_name','Calibrated Chlorophyll-Fluorescence Quality');
netcdf.putAtt(ncid,chlflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,chlflagid,'units','1');
netcdf.putAtt(ncid,chlflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,chlflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,chlflagid,'missing_value',-999);
netcdf.putAtt(ncid,chlflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,chlflagid,'flag_meanings','bad questionable good');

O2satid = netcdf.defVar(ncid,'O2_sat','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,O2satid,'long_name','Oxygen Saturation');
netcdf.putAtt(ncid,O2satid,'standard_name','fractional_saturation_of_oxygen_in_sea_water');
netcdf.putAtt(ncid,O2satid,'units','percent');
netcdf.putAtt(ncid,O2satid,'valid_range',[0 100]);
netcdf.putAtt(ncid,O2satid,'actual_range',[min(OXY_O2_sat(:)) max(OXY_O2_sat(:))]);
netcdf.putAtt(ncid,O2satid,'missing_value',-999);
netcdf.putAtt(ncid,O2satid,'ancillary_variables','O2_sat_qc');
netcdf.putAtt(ncid,O2satid,'origin','Computed');
netcdf.putAtt(ncid,O2satid,'notes','Computed using Gibbs SeaWater routine gsw_O2sol');

O2satflagid = netcdf.defVar(ncid,'O2_sat_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,O2satflagid,'long_name','Oxygen Saturation Quality');
netcdf.putAtt(ncid,O2satflagid,'standard_name','status_flag');
netcdf.putAtt(ncid,O2satflagid,'units','1');
netcdf.putAtt(ncid,O2satflagid,'valid_range',[-1 1]);
netcdf.putAtt(ncid,O2satflagid,'actual_range',[-1 1]);
netcdf.putAtt(ncid,O2satflagid,'missing_value',-999);
netcdf.putAtt(ncid,O2satflagid,'flag_values',[-1 0 1]);
netcdf.putAtt(ncid,O2satflagid,'flag_meanings','bad questionable good');

O2umolid = netcdf.defVar(ncid,'O2_umolkg','NC_DOUBLE',[depth_dimID pnum_dimID]);
netcdf.putAtt(ncid,O2umolid,'long_name','Oxygen concentration in micromolar per kilogram');
netcdf.putAtt(ncid,O2umolid,'standard_name','moles_of_oxygen_per_unit_mass_in_sea_water');
netcdf.putAtt(ncid,O2umolid,'units','micromole kilogram-1');
netcdf.putAtt(ncid,O2umolid,'valid_range',[0 max(OXY_O2_umolkg(:))]);
netcdf.putAtt(ncid,O2umolid,'actual_range',[min(OXY_O2_umolkg(:)) max(OXY_O2_umolkg(:))]);
netcdf.putAtt(ncid,O2umolid,'missing_value',-999);
netcdf.putAtt(ncid,O2umolid,'ancillary_variables','O2_umolkg_qc');
netcdf.putAtt(ncid,O2umolid,'origin','Computed');
netcdf.putAtt(ncid,O2umolid,'notes','Computed using Gibbs SeaWater routine gsw_O2sol');

O2umolflagid = netcdf.defVar(ncid,'O2_umolkg_qc','NC_DOUBLE',[depth_dimID pnum_dimID]);
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

% Profile number
netcdf.putVar(ncid,pnumid,1:size(CTD_time,2));

% Depth Grid
netcdf.putVar(ncid,depthid,1:size(CTD_time,1));

% Sensor ID
netcdf.putVar(ncid,sensorid,ID_sensor);

% Time
netcdf.putVar(ncid,tid,CTD_time);

% Longitude
netcdf.putVar(ncid,lonid,CTD_lon);

% Latitude
netcdf.putVar(ncid,latid,CTD_lat);

% in-situ Temperature
netcdf.putVar(ncid,Tid,CTD_T);

% in-situ Temperature flag
netcdf.putVar(ncid,Tflagid,floor(CTD_T_flag));

% Conservative temperature
netcdf.putVar(ncid,CTid,CTD_CT);

% Conservative temperature flag
netcdf.putVar(ncid,CTflagid,floor(CTD_CT_flag));

% Practical Salinity
netcdf.putVar(ncid,SPid,CTD_SP);

% Practical Salinity flag
netcdf.putVar(ncid,SPflagid,floor(CTD_SP_flag));

% Absolute Salinity
netcdf.putVar(ncid,SAid,CTD_SA);

% Absolute Salinity
netcdf.putVar(ncid,SAflagid,floor(CTD_SA_flag));

% In-situ density
netcdf.putVar(ncid,densid,CTD_rho);

% In-situ density
netcdf.putVar(ncid,densflagid,floor(CTD_rho_flag));

% bb470
netcdf.putVar(ncid,bb470id,FLS_bb470);

% bb470_flag
netcdf.putVar(ncid,bb470flagid,floor(FLS_bb470_flag));

% bb700
netcdf.putVar(ncid,bb700id,FLS_bb700);

% bb700_flag
netcdf.putVar(ncid,bb700flagid,floor(FLS_bb700_flag));

% chlorophyll
netcdf.putVar(ncid,chlid,FLS_chl);

% chlorophyll flag
netcdf.putVar(ncid,chlflagid,floor(FLS_chl_flag));

% Oxygen saturation
netcdf.putVar(ncid,O2satid,OXY_O2_sat);

% Oxygen saturation flag
netcdf.putVar(ncid,O2satflagid,floor(OXY_O2_sat_flag));

% Oxygen concentration
netcdf.putVar(ncid,O2umolid,OXY_O2_umolkg);

% Oxygen concentration flag
netcdf.putVar(ncid,O2umolflagid,floor(OXY_O2_umolkg_flag));

netcdf.close(ncid);
clear filename theind ncid *ID *id CTD*
