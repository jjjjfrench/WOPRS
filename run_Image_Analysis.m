function run_Image_Analysis(infilename, probetype, nChucks, projectname, threshold)
%% This function is the lead function in the Image_Analysis processing step.
%% 
%% Inputs: infilename, probetype, nChucks(max of 8), projectname, and threshold
%% (Threshold is only applicable to the CIP-G, and it represents the shading
%% threshold. It defaults to 50 if not provided.)
%%
%% Edited by Adam Majewski and Kevin Shaffer
%% 
%% Last edits made: 7/11/2019 by Kevin Shaffer
%%
%% Please consult the wiki for more detailed information and instructions

starpos = find(infilename == '*',1,'last');
slashpos = find(infilename == '/',1,'last');

if ~isempty(starpos)
    files = dir(infilename);
    filenums = length(files);
    filedir = infilename(1:slashpos);
else
    filenums = 1;
end

for i = 1:filenums
    if filenums > 1 || ~isempty(starpos)
        infilename = [filedir,files(i).name];
    end
    
    %  nChuck*nEvery should equal the total frame number 
    ncid = netcdf.open(infilename,'nowrite');
    time = netcdf.getVar(ncid, netcdf.inqVarID(ncid,'day'));
    nEvery = ceil(length(time)/nChucks);
    netcdf.close(ncid);
    numb=11:10+nChucks;  % Start from 11 to avoid sigle numbers in file name for later convinience

    % Assign the number of CPUs for this program
    if (nChucks > 1)
        parpool(nChucks)% Assign n CPUs to process
    end

    
    
    
    % Choose the start and end of chucks to be processed. Remember you can
    % split the chucks into different programs to process, since matlabpool can
    % only use 8 CPUs at once
    if (nChucks > 1)
        parfor iii=1:nChucks % 33:40  % iiith chuck will be processed 
            perpos = find(infilename == '.',1,'last');
            outfilename = [infilename(1:perpos-1),'_',num2str(iii),'.proc.cdf'];
            imgProc_sm(infilename,outfilename, probetype, iii, nEvery, projectname, threshold);  
        end
    else
        outfilename = [infilename(1:perpos-1),'.proc.cdf'];
        imgProc_sm(infilename,outfilename, probetype, iii, nEvery, projectname, threshold);
    end

    if (nChucks > 1)
        delete(gcp)
    end
end
end
