function [lambda,black,white] = model_fit()
%MODEL_FIT Summary of this function goes here
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

% Step 1: Determine location of files and files to process

% Step 1.1: Get list of files from selected folder
[pathname, filenames] = select_folder;

% Step 1.2: Check for background and calibration files, extract spectra if present
if cell2mat(strfind(filenames, 'black')) ~= 0 && cell2mat(strfind(filenames, 'white')) ~= 0
    [lambda, black, white] = get_bgcal(pathname, filenames); %% determine black & white spectra & lambda
else
    error('One or more of the calibration files is missing. (black.Master.Scope or white.Master.Scope)')
end

load('sbse_coeffs.mat')
load('model_params.mat')

% For each file,

for i=1:1 %length(filenames)
   current_filename = filenames{i};
   
   % Check if current file contains calibration or background measurement
   b = strfind(current_filename, 'black');
   w = strfind(current_filename, 'white');
   
   if isempty(b) && isempty(w) % If file doesn't contain cal or bg, continue
       % Step 2: Process the measured data in current file
       % Step 2.1: Retrieve the modified measured reflectance
       rmstar = calc_rmeas(pathname, current_filename, black, white);
       
       % Step 2.2: Correct for substitution error
       rm = corr_sbse(rmstar, sbse_coeffs);
       
       % Step 2.3: Determine uncertainty?
       
       % Step 3: Fit model to measured data
       %rdscat = calc_rd(conc,mua_param,musp)
       % Step 3.1: Determine coefficients using 'lsqcurvefit'
       % Step 3.2: Determine confidence interval using 'nlparci'
   else
       % Skip this file
   end
   

end

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

function [wl2, blk2, wht2] = get_bgcal(path, files)
% GET_BGCAL Extracts background and calibration spectra, and wavelength

% Retrieve black spec
fndblk = ones(1,length(files)) - double(cellfun('isempty',strfind(files, 'black')));
indblk = find(fndblk);
fileblk = strcat(path,files{indblk});
blkdata = dlmread(fileblk,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
intblk = dlmread(fileblk,' ', [6,3,6,3]); % reads in the integration time, space delimited
blk = (blkdata(:,2)/(intblk/1000))';
blk2=blk(453:1069);

% Retreive lambda values (500-705)
wl = blkdata(:,1)';
wl2=wl(453:1069);

% Retreive white spec
fndwht = ones(1,length(files)) - double(cellfun('isempty',strfind(files, 'white')));
indwht = find(fndwht);
filewht = strcat(path,files{indwht});
whtdata = dlmread(filewht,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
intwht = dlmread(filewht,' ', [6,3,6,3]); % reads in the integration time, space delimited
wht = (whtdata(:,2)/(intwht/1000))';
wht2=wht(453:1069);

end

function rmeas = calc_rmeas(path, file, bkg, cal)
% % CALC_RMEAS Retrieve the modified measured reflectance from Master.Scope
% % files
% 
curfile = strcat(path,file);
data = dlmread(curfile,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
inttime = dlmread(curfile,' ', [6,3,6,3]); % reads in the integration time, space delimited
spec = (data(:,2)/(inttime/1000))';
spec2 = spec(453:1069);
rmeas = (spec2 - bkg)./(cal - bkg);

end

function rcorr = corr_sbse(rmeas, coeff)

rcorr = ((coeff(:,1).*rmeas' + coeff(:,2))./(rmeas' + coeff(:,3)))';

end

function rdscat = calc_rd(beta,params,musp)
%RDSCAT Calculate diffuse reflectance with added scatter losses

% Step 1: Calculate diffuse reflectance using model
mua = beta(1)*beta(2)*params(1,:) + (beta(1)-beta(1)*beta(2))*params(2,:) + beta(3)*params(3,:) + beta(4)*params(4,:) + beta(5);

mutp = mua + musp;
mueff = sqrt(3.*mua.*mutp);
D = 1./(3.*mutp);

%A = 2.348; % for an nrel = 1.33 (for water)
A = 2.745; % for an nrel = 1.4 (for tissue)

rdcalc = musp./((mutp + mueff).*(1 + 2.*A.*D.*mueff));


end