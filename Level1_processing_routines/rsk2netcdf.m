clear
%%
% This code processes the raw data (Level 0 data) to produce individual
% netcdf files for each downcast and upcast, as well as for each
% instrument. Instruments include (1) the CTD, (2) the Oxygen sensor, and
% (3) the EcoPuck (backscatter, chlorophyll).

%%

% Set the counter for number of profiles recorded to 1
P_num = 1;

% List the RSK files located in the raw-data directory
path = './'
files = dir([path '*.rsk']);

% Looping through each RSK file in the folder.
for FF = 1:length(files)
    
    % Skip the whole 1st rsk file because the 2nd one includes the data from
    % the 1st one.
    if FF == 1
        continue
    end
    
    
    % Generate filename
    RSKfilename = fullfile(files(FF).folder,files(FF).name);
    
    % Open the RSK file using RBR Rusking toolbox
    RSK = RSKopen(RSKfilename);
    % Extract the data from RSK file
    RSK = RSKreaddata(RSK);
    
    
    % Record Section number from filename
    disp(['Section number: ',RSKfilename(end-5:end-4)]);
    Sec_num = str2double(RSKfilename(end-5:end-4));
    
    % Correct for the issue if trying to read a 1-digit number (Sec_num<10)
    if isnan(Sec_num)==1
        Sec_num = str2double(RSKfilename(end-4));
    end
    if isnan(Sec_num)==1
        error('Could not find Section number')
    end
    
    %% PLOT PROFILES
%     figure
%     scatter(RSK.data.tstamp,...
%         RSK.data.values(:,3),...
%         35,'ok','filled')
%     hold on
%     
%     for ii = [1:length(RSK.profiles.downcast.tstart)]
%         scatter(RSK.data.tstamp(RSK.data.tstamp>=RSK.profiles.downcast.tstart(ii) & RSK.data.tstamp<=RSK.profiles.downcast.tend(ii)),...
%             RSK.data.values(RSK.data.tstamp>=RSK.profiles.downcast.tstart(ii) & RSK.data.tstamp<=RSK.profiles.downcast.tend(ii),3),...
%             35,'.r')
%         hold on
%     end; clear ii
%     
%     for ii = [1:length(RSK.profiles.upcast.tstart)]
%         scatter(RSK.data.tstamp(RSK.data.tstamp>=RSK.profiles.upcast.tstart(ii) & RSK.data.tstamp<=RSK.profiles.upcast.tend(ii)),...
%             RSK.data.values(RSK.data.tstamp>=RSK.profiles.upcast.tstart(ii) & RSK.data.tstamp<=RSK.profiles.upcast.tend(ii),3),...
%             35,'.g')
%         hold on
%     end; clear ii
%     datetick('x','ddmmm HH:MM')
%     set(gca,'fontsize',16)
%     ylabel('Depth [m]')
%     grid on
    
    %% Generate NetCDF Profile file
    
    % Loops through profiles.
    for ii = 1:length(RSK.profiles.downcast.tstart)
        
        % Based on manual vizualization, some profiles were rejected. Once
        % the indices were established, the process was automated by
        % removing those profiles automatically.
        
        % Find data indices for the considered profile, for both upcasts
        % and downcasts.
        ind_down = find(RSK.data.tstamp>=RSK.profiles.downcast.tstart(ii) &...
            RSK.data.tstamp<=RSK.profiles.downcast.tend(ii));
        ind_up = find(RSK.data.tstamp>=RSK.profiles.upcast.tstart(ii) &...
            RSK.data.tstamp<=RSK.profiles.upcast.tend(ii));
        
        %% Plots the profile for approval
        
        figure
        % Section number 16 includes a known offset on the the timestamp of
        % the data. Offset was identified during the cruise and confidently
        % quantified.
        if Sec_num == 16
            scatter(RSK.data.tstamp(ind_down)+6724.46374999,RSK.data.values(ind_down,3),35,'.r')
            hold on
            scatter(RSK.data.tstamp(ind_up)+6724.46374999,RSK.data.values(ind_up,3),35,'.g')
        else
            scatter(RSK.data.tstamp(ind_down),RSK.data.values(ind_down,3),35,'.r')
            hold on
            scatter(RSK.data.tstamp(ind_up),RSK.data.values(ind_up,3),35,'.g')
        end
        ylabel('Depth [m]'); datetick('x','ddmmm HH:MM'); legend('downcast','upcast','location','best')
        grid on; set(gca,'fontsize',16);axis ij
        
        %Ask user if profile should be printing into NetCDF file
        %         answer = questdlg(['Should the cast #',num2str(P_num),' in Section #', num2str(Sec_num),' be recorded in NetCDF file? (ii = ',num2str(ii),')'], ...
        %            'Cast approval', ...
        %            'Yes','No','Cancel','Yes');
        answer = 'Yes';
        
        % Handle response
        switch answer
            % No file is produce
            case 'No'
                close
                continue
                
                % abort the code
            case 'Cancel'
                error('Aborted by user')
                
                % Print data into file
            case 'Yes'
                % close the figure
                close
                
                %% CTD
                
                % Create the NetCDF file for CTD data
                disp('Saving CTD data ...')
                %run ctd_netcdf
                
                % Create the NetCDF file for Oxygen data
                disp('Saving Oxygen data ...')
                %run oxy_netcdf
                
                % Create the NetCDF file for Optical data
                disp('Saving Optical data ...')
                run fls_netcdf
                
                P_num = P_num+1;
        end
        
    end
end
