function Generate_SizeDist(PROC_files,outfile,tas,timehhmmss,probename,num_rejects,IA_threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is called by run_SizeDist to generate NetCDF files
% containing size distribution data. 
%
% The NetCDF files produced by this function include the following
% information:
% -concentration of particles by rejection type, seperated into bins based
% on diameter size
% -concentrations of particles separated into bins based on roundness
% and aspect ratio values
% -particle counts for accepted particles seperated into bins based on
% diameter size
% -particle counts for rejected particles seperated into bins based on
% diameter size
% -rejection ratio across all bins for each probe
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    PROC_files struct
    outfile (1,:) string
    tas (1,:) double
    timehhmmss (1,:) double
    probename (1,:) char
    num_rejects (1,1) {mustBeNumeric,mustBeReal}
    IA_threshold (1,1) {mustBeNumeric,mustBeReal}=10; %Default interarrival threshold (clock cycles) used if no changes are made to run_SizeDist and Generate_SizeDist is called w/o arg
end 
%% Begin generating the size distribution
disp(['Generating Size Distribution for the ',probename]);

%% Read the setup file with bin edges
[roundness_bin_edges,num_round_bins,In_status,in_status_value,diam_bin_edges,num_diam_bins,mid_bin_diams,aspect_ratio_bin_edges,num_ar_bins,circularity_bin_edges,num_circ_bins]=setup_SizeDist(probename);

%% Create NaN arrays
accepted_counts  = NaN(length(timehhmmss),num_diam_bins);
total_accepted_counts = NaN(length(timehhmmss),1);
total_rejected_counts = NaN(length(timehhmmss),num_rejects+1);
roundness_counts = NaN(num_diam_bins,num_round_bins,length(timehhmmss));
circularity_counts = NaN(num_diam_bins,num_circ_bins,length(timehhmmss));
aspect_ratio_counts = NaN(num_diam_bins,num_ar_bins,length(timehhmmss));
sample_volume = NaN(length(timehhmmss),num_diam_bins);
size_dist = NaN(length(timehhmmss),num_diam_bins);

switch probename
    case '2DS'
        accepted_counts_H = NaN(length(timehhmmss),num_diam_bins);
        accepted_counts_V = NaN(length(timehhmmss),num_diam_bins);
        size_dist_H = NaN(length(timehhmmss),num_diam_bins);
        size_dist_V = NaN(length(timehhmmss),num_diam_bins);
end
min_bin_diams = diam_bin_edges(1:end-1);
max_bin_diams = diam_bin_edges(2:end);

%***************************************************************************************************************
%***************************************************************************************************************

%% Loop over each PROC file in the directory
for x = 1:length(PROC_files)
    
    % Define input file and initialize time variable
    disp(['Reading in file: ',PROC_files(x).folder,'/',PROC_files(x).name])
    infile = netcdf.open([PROC_files(x).folder,'/',PROC_files(x).name],'nowrite');
    times = netcdf.getVar(infile,netcdf.inqVarID(infile,'Time'));

    % Read in other variables
    diameter = netcdf.getVar(infile,netcdf.inqVarID(infile,'diameter'));
    artifacts = netcdf.getVar(infile,netcdf.inqVarID(infile,'artifact_status'));
    roundness = netcdf.getVar(infile,netcdf.inqVarID(infile,'roundness'));
    circularity = netcdf.getVar(infile,netcdf.inqVarID(infile,'circularity'));
    aspect_ratio = netcdf.getVar(infile,netcdf.inqVarID(infile,'aspect_ratio'));
    inter_arrival = netcdf.getVar(infile,netcdf.inqVarID(infile,'inter_arrival'));
    shattered = inter_arrival <= IA_threshold; %logical indices of all shattered particles
    ind_accepted = (artifacts==1); %logical indices of accepted (1) artifacts flag
    %for accepted particles, overwrites artifact status to new "shattered"
    %artifact of max(artifacts)+1, otherwise status remains accepted (1):
    artifacts(ind_accepted) = (shattered(ind_accepted)*max(artifacts))+1;
    switch probename
        case '2DS'
            channel = netcdf.getVar(infile,netcdf.inqVarID(infile,'channel'));
    end
    in_status = netcdf.getVar(infile,netcdf.inqVarID(infile,'in_status'));
    
    % Read in global attributes necessary for sample area calculation
    diodesize = netcdf.getAtt(infile, netcdf.getConstant('NC_GLOBAL'),'Diode size');
    num_diodes = netcdf.getAtt(infile, netcdf.getConstant('NC_GLOBAL'),'Number of diodes');
    armdst = netcdf.getAtt(infile, netcdf.getConstant('NC_GLOBAL'),'Arm distance');
    wavelength = netcdf.getAtt(infile, netcdf.getConstant('NC_GLOBAL'),'Wavelength');
    num_rejects = netcdf.getAtt(infile, netcdf.inqVarID(infile,'artifact_status'),'Number of artifact statuses')+1;
    
    %% Create the output file
    [f,varid]=define_outfile_SizeDist(probename,num_rejects,timehhmmss,outfile,num_round_bins,num_diam_bins,In_status,num_ar_bins,num_circ_bins,IA_threshold);
    
    % Times shared by aircraft and probe data
    times_both = times(ismember(times,timehhmmss));
    
    % Increment PROC if no matching times
    if isempty(times), continue, end
    
    first_PROC_time = times_both(1);
    last_PROC_time = times_both(end);

    % Find the first and last time that is available in both the PROC and
    % aircraft files. We will loop over these times to create the SD file.
    first_time_index = find(timehhmmss == first_PROC_time);
    end_time_index = find(timehhmmss == last_PROC_time);

    %% Calculate sample area depending on whether we are considering center-in or only all-in images 
    
    % Calculate for each bin
    switch In_status
        case {'Center-in','center-in','Centerin','centerin','Center','center'}
            for i=1:length(mid_bin_diams)
                DOF = (8*mid_bin_diams(i)^2) / (4*wavelength);%Calculate depth-of-field
                along_beam = min(armdst,DOF);%If the DOF is larger than the distance between the probe arms, then the DOF is equal to the distance between probe arms
                sample_area(i) = (along_beam * diodesize * (num_diodes-2+1)) / 100; %Convert to cm^2
            end
        case {'All-in','all-in','Allin','allin','All','all'}
            for i=1:length(mid_bin_diams)
                DOF = (8*mid_bin_diams(i)^2) / (4*wavelength);%Calculate depth-of-field
                along_beam = min(armdst,DOF);%If the DOF is larger than the distance between the probe arms, then the DOF is equal to the distance between probe arms
                sample_area(i) = (along_beam * diodesize * (num_diodes-2+1) - mid_bin_diams(i)) / 100; %Convert to cm^2
            end
    end



    %% Loop over every second
    for i = first_time_index:end_time_index
        current_hhmmss = timehhmmss(i);

        % Find the indices of all the particles that correspond to this time
        good_indices = find(times == current_hhmmss);
        
        % Set the counts to 0 at every time when there is both an aircraft
        % file and a PROC file. NaN's otherwise
        accepted_counts(i,:) = 0;
        total_accepted_counts(i) = 0;
        total_rejected_counts(i,:) = 0;
        roundness_counts(:,:,i) = 0;
        circularity_counts(:,:,i) = 0;
        aspect_ratio_counts(:,:,i) = 0;
        switch probename
            case '2DS'
                accepted_counts_H(i,:) = 0;
                accepted_counts_V(i,:) = 0;
                size_dist_H = zeros(length(timehhmmss),num_diam_bins)*NaN;
                size_dist_V = zeros(length(timehhmmss),num_diam_bins)*NaN;
        end
        size_dist(i,:) = 0;
        
        %% Calculate sample volume second-by-second using tas from aircraft data
        sample_volume(i,:) = (tas(i)*100) * sample_area; %In cm^3
        
        % Read in the area ratio and artifact status for
        % the particles at this time. Then add to appropriate counts
        if ~isempty(good_indices) % If there is more than 1 particle in this second, continue
            good_diameters = diameter(good_indices); % Find the diameters of all particles in this time
            good_artifacts = artifacts(good_indices); % Find the artifact statuses
            good_roundness = roundness(good_indices); % Find the area ratios
            good_circularity = circularity(good_indices); % Find the circularities
            good_aspect_ratio = aspect_ratio(good_indices); %Find the aspect ratios
            switch probename
                case '2DS'
                    good_channel = channel(good_indices); % Find the H or V status
            end
            good_in_status = in_status(good_indices); % Find if the image is center-in, center-out, or all-in
            
            for j=1:length(good_diameters)
                %% Rejection Counts
                % Find which bin the given diameter fits into
                switch good_in_status(j)
                    case in_status_value
                        switch good_artifacts(j)
                            case 1 % Particle is not rejected
                                %% Rejection counts
                                lower_bin_numbers = find(diam_bin_edges <= good_diameters(j));
                                diam_bin_number = lower_bin_numbers(end); % This is the number of the bin that the diameter fits into
                                if diam_bin_number == length(diam_bin_edges) continue, end % Advances particle if D doesnt fit into a bin
                                accepted_counts(i,diam_bin_number) = accepted_counts(i,diam_bin_number) + 1;
                                total_accepted_counts(i) = total_accepted_counts(i) + 1;
                                switch probename
                                    case '2DS'
                                        switch good_channel(j)
                                            case 'H'
                                                accepted_counts_H(i,diam_bin_number) = accepted_counts_H(i,diam_bin_number) + 1;
                                            case 'V'
                                                accepted_counts_V(i,diam_bin_number) = accepted_counts_V(i,diam_bin_number) + 1;
                                            otherwise
%                                                 disp('Error')
%                                                 disp('Shutting down')
                                        end
                                end
                        
                                %% Roundness counts
                                if ~isnan(good_roundness(j))
                                    lower_bin_numbers = find(roundness_bin_edges <= good_roundness(j));
                                    round_bin_number = lower_bin_numbers(end); % This is the number of the bin that the roundness fits into
                                    % Add counts to roundness bins
                                    roundness_counts(diam_bin_number,round_bin_number,i) = roundness_counts(diam_bin_number,round_bin_number,i) + 1;
                                end
                                
                                if ~isnan(good_circularity(j))
                                    %% Circularity counts
                                    lower_bin_numbers = find(circularity_bin_edges <= good_circularity(j));
                                    circ_bin_number = lower_bin_numbers(end); % This is the number of the bin that the circularity fits into
                                    % Add counts to roundness bins
                                    circularity_counts(diam_bin_number,circ_bin_number,i) = circularity_counts(diam_bin_number,circ_bin_number,i) + 1;
                                end

                                if ~isnan(good_aspect_ratio(j))
                                    %% Aspect ratio counts
                                    lower_bin_numbers = find(aspect_ratio_bin_edges <= good_aspect_ratio(j));
                                    ar_bin_number = lower_bin_numbers(end); % This is the number of the bin that the area ratio fits into
                                    % Add counts to roundness bins
                                    aspect_ratio_counts(diam_bin_number,ar_bin_number,i) = aspect_ratio_counts(diam_bin_number,ar_bin_number,i) + 1;
                                end
                                
                            otherwise
                                total_rejected_counts(i,good_artifacts(j)) = total_rejected_counts(i,good_artifacts(j)) + 1;
                        end             
                end
            end
        end
    end
end

% Calculate size distributions
size_dist = accepted_counts ./ sample_volume;
switch probename
    case '2DS'
        size_dist = size_dist / 2; %Need to divide concentrations by 2 becuase the 2DS has two channels
        size_dist_H = accepted_counts_H ./ sample_volume;
        size_dist_V = accepted_counts_V ./ sample_volume;
end

% Normalize by bin width
for i=1:num_diam_bins
    size_dist(:,i) = accepted_counts(:,i) ./ (max_bin_diams(i) - min_bin_diams(i));
    switch probename
        case '2DS'
            size_dist_H(:,i) = accepted_counts_H(:,i) ./ (max_bin_diams(i) - min_bin_diams(i));
            size_dist_V(:,i) = accepted_counts_V(:,i) ./ (max_bin_diams(i) - min_bin_diams(i));
    end
end
        
% Assign data to outfile variables
netcdf.putVar ( f, varid.time, timehhmmss);
netcdf.putVar ( f, varid.Accepted_counts, accepted_counts);
netcdf.putVar ( f, varid.total_accepted_counts, total_accepted_counts);
netcdf.putVar ( f, varid.bin_min, min_bin_diams);
netcdf.putVar ( f, varid.bin_max, max_bin_diams);
netcdf.putVar ( f, varid.bin_mid, mid_bin_diams);
netcdf.putVar ( f, varid.roundness_counts, roundness_counts);
netcdf.putVar ( f, varid.circularity_counts, circularity_counts);
netcdf.putVar ( f, varid.aspect_ratio_counts, aspect_ratio_counts);
netcdf.putVar ( f, varid.total_reject_counts, total_rejected_counts);       
switch probename
    case '2DS'
        netcdf.putVar ( f, varid.Accepted_counts_H, accepted_counts_H);
        netcdf.putVar ( f, varid.Accepted_counts_V, accepted_counts_V);
        netcdf.putVar ( f, varid.size_dist_2DS_H, size_dist_H);
        netcdf.putVar ( f, varid.size_dist_2DS_V, size_dist_V);
end
netcdf.putVar ( f, varid.size_dist, size_dist);
netcdf.putVar ( f, varid.sample_volume, sample_volume);
netcdf.close(f);

end