function [roundness_bin_edges,num_round_bins,In_status,in_status_value,diam_bin_edges,num_diam_bins,mid_bin_diams,aspect_ratio_bin_edges,num_aspect_ratio_bins,circularity_bin_edges,num_circ_bins]=setup_SizeDist(probename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function sets all the bin edges for generating size distributions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% All probetypes:
roundness_bin_edges = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
aspect_ratio_bin_edges = [1 1.2 1.4 1.6 1.8 2 2.25 2.5 2.75 3 4 5 7.5 10 100 Inf];
circularity_bin_edges = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 Inf];

%The following should have a value of ‘center-in’ or ‘all-in’
In_status = 'all-in';

%% Probetype-specific:
switch probename
    case '2DC'
        res = 25; %um pixel resolution
        n = 64;
    case '2DP'
        res = 200; %um pixel resolution
        n = 32;
    case 'HVPS'
        res = 150; %um pixel resolution
        n = 128;
    case '2DS'
        res = 10; %um pixel resolution
        n = 128;
    case 'CIP'
        res = 25; %um pixel resolution
        n = 64;
    otherwise 
        disp('ERROR: Probetype is not supported. Please enter one of the following: 2DP, 2DC, 2DS, HVPS, or CIP. Note: Matlab is case sensitive')
        return;
end
diam_bin_edges = linspace(res,(2*n)*res,2*n); % generates bins of up to 2 times the width of the photodiode array (for center out)

switch In_status
    case {'Center-in','center-in','Centerin','centerin','Center','center'}
        in_status_value = {'A','I'}; %All-in or Center-in
    case {'All-in','all-in','Allin','allin','All','all'}
        in_status_value = {'A'}; %All-in only
    otherwise 
        disp('ERROR: In_status choice does not make sense. Check to make sure center-in or all-in is chosen in the setup file.')
        return;
end

num_round_bins = length(roundness_bin_edges) -1;
num_aspect_ratio_bins = length(aspect_ratio_bin_edges) -1;
num_diam_bins = length(diam_bin_edges) -1;
num_circ_bins = length(circularity_bin_edges)-1;

for i=1:num_diam_bins
    mid_bin_diams(i) = mean(diam_bin_edges(i:i+1))/1000; %Get mid_bin_diams in millimeters for sample area calculation
end

end