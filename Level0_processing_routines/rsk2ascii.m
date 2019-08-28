clear

%% USER DEFINED PARAMTERS
% Path to RSK-files
rsk_path = '../RSK_files';

% Path to save ASCII-file
ascii_path = pwd;

% Add the folder with the necessary dependencies to the path
addpath([ascii_path,'/dependencies/']);

% Personalized header
persoheader = 1;    % Set to 1 if header is prescribed by user
myheader = {'# UTC Timestamp [milliseconds since January 1 1970]',...
    'Conductivity [mS/cm]',...
    'Temperature [Degrees_C]',...
    'Pressure [dbar]',...
    'Dissolved O2 saturation [percent]',...
    'Backscatter at 470nm [counts]',...
    'Backscatter at 700nm [counts]',...
    'Fluorescence [counts]'};

%% CORE CODE

% List rsk files in folder
files = dir([rsk_path,'/*.rsk']);

% Loop through the RKS file and print them as ASCII
for ii = 22%1:length(files)
    
    databasefile = fullfile(files(ii).folder,files(ii).name);

    % opens the database
    mksqlite('open', databasefile);
    
    % Get variable names
    instru = mksqlite('SELECT * from instruments');    
    channels = mksqlite('SELECT * from channels');
    data = mksqlite('SELECT * from data');
    
    % Create the ASCII file
    fid = fopen([ascii_path,'/EcoCTD_S',num2str(ii),'.ascii'],'w');
    
    %% HEADER
    % Prints info about the RBR data logger
    fprintf(fid,'%s\n','# BEGIN HEADER');
    fprintf(fid,'%s\n','# Header has 14 lines');    
    fprintf(fid,'%s\n','# This file follows the guidelines for ASCII-files format in Earth Sciences: earthdata.nasa.gov/user-resources/standards-and-references/ascii-file-format-guidelines-for-earth-science-data');    
    fprintf(fid,'%s\n','# Cruise Info: CALYPSO (ONR) pilot cruise 25/05/2018-03/06/2018');    
    fprintf(fid,'%s\n','# Principal Investigator: Dr. A. Mahadevan (amahadevan@whoi.edu)');    
    fprintf(fid,'%s\n',['# ASCII file generated from ',files(ii).name]);    
    fprintf(fid,'%s\n',['# ASCII file generated on ',date]);    
    fprintf(fid,'%s\n',['# RBR Instrument model: ',instru.model]);
    fprintf(fid,'%s%d\n','# RBR Instrument Serial Number: ',instru.serialID);
    fprintf(fid,'%s\n',['# RBR Instrument Firmware version: ',instru.firmwareVersion]);
    if ii == 16
        fprintf(fid,'%s\n',['# Data recorded from ',...
        datestr(data(1).tstamp/1000/86400+datenum([1970 1 1 0 0 0])+6724.46374999),...
        ' to ',...
        datestr(data(end).tstamp/1000/86400+datenum([1970 1 1 0 0 0])+6724.46374999)]);
    else
        fprintf(fid,'%s\n',['# Data recorded from ',...
        datestr(data(1).tstamp/1000/86400+datenum([1970 1 1 0 0 0])),...
        ' to ',...
        datestr(data(end).tstamp/1000/86400+datenum([1970 1 1 0 0 0]))]);
    end
    fprintf(fid,'%s\n','# Missing Value = NaN');

    
    % Record channel number that is directly measured
    measuredchannelnum = [];
    for jj = 1:length(channels)
        % only print if quantity is directly measured
        if channels(jj).isMeasured==1
            measuredchannelnum = cat(1,measuredchannelnum,jj);
        end
    end
    
    % loop through the directly measured channels
    if persoheader == 1
        for jj = 1:length(myheader)
            if jj == length(myheader)
                fprintf(fid,'%s\n',myheader{jj});
            else
                fprintf(fid,'%s, ',myheader{jj});
            end
        end; clear jj
    else
        for jj = 0:max(measuredchannelnum)
            if jj == 0
                fprintf(fid,'%s, ','# UTC Timestamp [milliseconds since January 1 1970]');
            elseif jj == max(measuredchannelnum)
                fprintf(fid,'%s\n',[channels(jj).longNamePlainText,' [',channels(jj).unitsPlainText,']']);
            else
                fprintf(fid,'%s, ',[channels(jj).longNamePlainText,' [',channels(jj).unitsPlainText,']']);
            end
        end; clear jj
    end
    fprintf(fid,'%s\n','# END HEADER');

    %% DATA
    % Get the data from database
    tic
    
    % Converts from database to matrix
    % Get fieldnames in data structure
    names = fieldnames(data);
    % Create empty matrix
    thedata = NaN*zeros(length(data),length(measuredchannelnum)+1);
    
    % extract the timestamp
    for jj = 1:length(data)
        if isempty(eval('data(jj).tstamp;'))
            continue
        else
            eval(['thedata(jj,1) = data(jj).',names{1},';']);
        end
    end; clear jj
    
    if ii == 16
        thedata(:,1) = thedata(:,1)+6724.46374999*86400*1000;
    end
    % extract the directly measured variables only
    for ff = 1:length(measuredchannelnum)
        for jj = 1:length(data)
            if isempty(eval(['data(jj).',names{measuredchannelnum(ff)+1},';']))
                %warning('empty val')
                continue
            else
                eval(['thedata(jj,measuredchannelnum(ff)+1) = data(jj).',names{measuredchannelnum(ff)+1},';']);
            end
        end; clear jj
    end; clear ff
    toc
    
    % Print into the ASCII file
    fprintf(fid,'%f, %f, %f, %f, %f, %d, %d, %d\n',thedata');
    % close ascii file
    fclose(fid)
    
    disp('File written !!')
    
    % Close the database
    mksqlite('close')    
end