function run_SizeDist(ncfile,PROC_directory,probe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses the nc file for a given flight and uses the corresponding PROC files
% for a given probe to generate second-by-second size distributions.
%
% Example Inputs:
% ncfile - '/kingair_data/snowie17/work/123456.c1.nc'
% PROC_directory - '/kingair_data/snowie17/OAP_processed/123456'
% probe - '2DS'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_nc = netcdf.open(ncfile);
timesec = netcdf.getVar(main_nc, netcdf.inqVarID(main_nc, 'time'), 'double');
try 
    tas = netcdf.getVar(main_nc, netcdf.inqVarID(main_nc, 'tas'),'double');
catch 
    tas = netcdf.getVar(main_nc, netcdf.inqVarID(main_nc, 'TAS'),'double');
end
netcdf.close(main_nc);
days = floor(timesec/86400);
hours = floor(mod(timesec,86400)/3600);
minutes = floor(mod(timesec,3600)/60);
seconds = mod(timesec,60);
timehhmmss = mod(seconds+minutes*100+hours*10000,240000); % Rolls to 0 at 24:00:00UTC after start day

% Find the date and project directory from the nc file
c1_name = split(ncfile,'/');
c1_name = char(c1_name(end));
parts = split(c1_name,{'.','_','-'}); %split nc filename string on '.', '_' or '-' delimiters
nums = regexp(parts,'\d+','match'); %only numerals
allnums = vertcat(nums{:}); %fix nested cells
len = cellfun(@length,allnums); %length of nums cell strings
date = char(allnums(len == max(len))); %date set to continuous set of numerals of max length
PROC_directory = ([PROC_directory,'/']);

% In-status to append to name
[~,~,in_status,in_status_value,~,~,~,~,~]=setup_SizeDist(probe);

% Sanity check to make sure we have true airspeed data for every second
if length(timehhmmss)==length(tas)
    disp(['Length of the nc file is: ',num2str(length(tas)),' seconds'])
else 
    disp('ERROR: Length of the time and true airspeed variables is not equal. Shutting down.')
    return;
end

% Find the appropriate file(s) based on the probetype and then run
% Generate_SizeDist.
switch probe
    case {'2DS','2ds'}

        files = dir([PROC_directory,'PROC.*.2DS.cdf']);
        filenums = length(files);
        if filenums < 1
            disp('No 2DS files found in the provided directory')
            return;
        end
        
        % Read in the artifact statuses to see how many there are. This
        % needs to be done so that we can set one of the dimension sizes in
        % Generate_SizeDist. By reading this in here, artifact statuses can
        % be added, without affecting the size distribution code.
        PROC_file = netcdf.open([PROC_directory,files(1).name],'nowrite');
        num_rejects = netcdf.getAtt(PROC_file, netcdf.inqVarID(PROC_file,'artifact_status'),'Number of artifact statuses');
        IA_threshold = 1e-5;
        inFile = [PROC_directory,files(1).name];
        proc_pos = find(inFile == 'P',1,'last');
        outFile = [inFile(1:proc_pos-1),'SD.',date,'_',upper(in_status(1)),'.2DS.cdf']
        Generate_SizeDist(files,outFile,tas,floor(timehhmmss),'2DS',num_rejects,IA_threshold);    
    case {'HVPS','hvps'}
        files = dir([PROC_directory,'PROC.*.HVPS.cdf']);
        filenums = length(files);
        if filenums < 1
            disp('No HVPS files found in the provided directory')
            return;
        end
        
        % Read in the artifact statuses to see how many there are. This
        % needs to be done so that we can set one of the dimension sizes in
        % Generate_SizeDist. By reading this in here, artifact statuses can
        % be added, without affecting the size distribution code.
        PROC_file = netcdf.open([PROC_directory,files(1).name],'nowrite');
        num_rejects = netcdf.getAtt(PROC_file, netcdf.inqVarID(PROC_file,'artifact_status'),'Number of artifact statuses');
        IA_threshold = 1e-4;
        inFile = [PROC_directory,files(1).name];
        proc_pos = find(inFile == 'P');
        if length(proc_pos) > 1
            proc_pos = proc_pos(end-1);
        else
            proc_pos = proc_pos(end);
        end
        outFile = [inFile(1:proc_pos-1),'SD.',date,'_',upper(in_status(1)),'.HVPS.cdf']
        Generate_SizeDist(files,outFile,tas,floor(timehhmmss),'HVPS',num_rejects,IA_threshold);
        
    case {'CIPG','cip','cipg','CIP'}

        files = dir([PROC_directory,'PROC.*.CIP.cdf']);
        filenums = length(files);
        if filenums < 1
            disp('No CIP files found in the provided directory')
            return;
        end
        
        % Read in the artifact statuses to see how many there are. This
        % needs to be done so that we can set one of the dimension sizes in
        % Generate_SizeDist. By reading this in here, artifact statuses can
        % be added, without affecting the size distribution code.
        PROC_file = netcdf.open([PROC_directory,files(1).name],'nowrite');
        num_rejects = netcdf.getAtt(PROC_file, netcdf.inqVarID(PROC_file,'artifact_status'),'Number of artifact statuses');
        IA_threshold = 1e-5;
        inFile = [PROC_directory,files(1).name];
        proc_pos = find(inFile == 'P');
        if length(proc_pos) > 1
            proc_pos = proc_pos(end-1);
        else
            proc_pos = proc_pos(end);
        end
        outFile = [inFile(1:proc_pos-1),'SD.',date,'_',upper(in_status(1)),'.CIP.cdf']
        Generate_SizeDist(files,outFile,tas,floor(timehhmmss),'CIP',num_rejects,IA_threshold);

    case {'2DP','2dp'}

        files = dir([PROC_directory,'PROC.*.2DP.cdf']);
        filenums = length(files);
        if filenums < 1
            disp('No 2DP files found in the provided directory')
            return;
        end
        
        % Read in the artifact statuses to see how many there are. This
        % needs to be done so that we can set one of the dimension sizes in
        % Generate_SizeDist. By reading this in here, artifact statuses can
        % be added, without affecting the size distribution code.
        PROC_file = netcdf.open([PROC_directory,files(1).name],'nowrite');
        num_rejects = netcdf.getAtt(PROC_file, netcdf.inqVarID(PROC_file,'artifact_status'),'Number of artifact statuses');
        IA_threshold = 1e-3;
        inFile = [PROC_directory,files(1).name];
        proc_pos = find(inFile == 'P');
        if length(proc_pos) > 1
            proc_pos = proc_pos(end-1);
        else
            proc_pos = proc_pos(end);
        end
        outFile = [inFile(1:proc_pos-1),'SD.',date,'_',upper(in_status(1)),'.2DP.cdf']
        Generate_SizeDist(files,outFile,tas,floor(timehhmmss),'2DP',num_rejects,IA_threshold);
        
    case {'2DC','2dc'}
        
        files = dir([PROC_directory,'PROC.*.2DC.cdf']);
        filenums = length(files);
        if filenums < 1
            disp('No 2DC files found in the provided directory')
            return;
        end
        
        % Read in the artifact statuses to see how many there are. This
        % needs to be done so that we can set one of the dimension sizes in
        % Generate_SizeDist. By reading this in here, artifact statuses can
        % be added, without affecting the size distribution code.
        PROC_file = netcdf.open([PROC_directory,files(1).name],'nowrite');
        num_rejects = netcdf.getAtt(PROC_file, netcdf.inqVarID(PROC_file,'artifact_status'),'Number of artifact statuses');
        
        inFile = [PROC_directory,files(1).name];
        proc_pos = find(inFile == 'P',1,'last');
        outFile = [inFile(1:proc_pos-1),'SD.',date,'_',upper(in_status(1)),'.2DC.cdf']
        Generate_SizeDist(files,outFile,tas,floor(timehhmmss),'2DC',num_rejects);
        
end

end
