function [] = model_fit()
%MODEL_FIT Summary of this function goes here
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, filenames] = select_folder;

% Step 1.2: Check for background and calibration files
if exist(cell2mat(strfind(filenames, 'black')))~=0 && exist(cell2mat(strfind(filenames, 'white')))~=0
    %determine black & white spectra
else
    error('One or more of the calibration files is missing. (black.Master.Scope or white.Master.Scope)')
end
% Step 1.3: Get lambda
load('wavelength.mat')

% For each file,
i=1;
%for i=1:length(filenames)
    current_filename = filenames{i};
    
    % Check if current file contains calibration or background measurement
    b = strfind(current_filename, 'black');
    w = strfind(current_filename, 'white');
    
    if isempty(b) && isempty(w) % If file doesn't contain cal or bg, continue
        % Step 2: Process the measured data in current file
        % Step 2.1: Retrieve the modified measured reflectance
        rmstar = calc_rmeas(pathname, current_filename);
        
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

function rmeas = calc_rmeas(path, file)
% CALC_RMEAS Retrieve the modified measured reflectance from Master.Scope
% files

file = strcat(pathname,filename);
data = dlmread(file,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
int_time = dlmread(file,' ', [6,3,6,3]); % reads in the integration time, space delimited
spec = (data(:,2)/(int_time/1000))';

end