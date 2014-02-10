function [] = model_fit()
%MODEL_FIT Summary of this function goes here
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, filenames] = select_folder;

% For each file,
i=1;
%for i=1:length(filenames)
    current_file=strcat(pathname,filenames{i});
    
    % check if current file contains calibration or background measurement
    b = strfind(current_file, 'black');
    w = strfind(current_file, 'white');
    
    if isempty(b) && isempty(w) % if file doesn't contain cal or bg, continue
        % Step 2: Process the measured data in current file
        % Step 2.1: Retrieve the modified measured reflectance
        % Step 2.2: Correct for substitution error
        % Step 2.3: Determine uncertainty?
        
        % Step 3: Fit model to measured data
        % Step 3.1: Determine coefficients using 'lsqcurvefit'
        % Step 3.2: Determine confidence interval using 'nlparci'
    else
        % Skip this file
    end
%end




end

function [pathcat, strucell2] = select_folder
%SELECT_FOLDER Get a list of files from selected folder

folder_name = uigetdir;
pathcat = strcat(folder_name,'\');
files = dir(folder_name);

strucell = struct2cell(files);
strucell1 = strucell(1,:);
strucell2 = strucell1(3:numel(files));

end