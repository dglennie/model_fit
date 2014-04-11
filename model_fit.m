function [param, eicorr] = model_fit()
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

m = 1;

% Initial set up for fitting
load('sbse_coeffs.mat')
load('model_params.mat')
param0 = [0.0076, 0.8421, -0.0017, 5.0798, 0.0321];
lb = [0, 0.5, -0.01, 0, -0.2];
ub = [0.02, 1, 0.01, 100, 0.2];
options = optimset('Algorithm','levenberg-marquardt');

% For each file,

for i=4:4 %length(filenames)
   current_filename = filenames{i};
   xlsfilenames{m} = current_filename(1:end-13);
   
   % Check if current file contains calibration or background measurement
   b = strfind(current_filename, 'black');
   w = strfind(current_filename, 'white');
   
   if isempty(b) && isempty(w) % If file doesn't contain cal or bg, continue
       % Step 2: Process the measured data in current file
       % Step 2.1: Retrieve the modified measured reflectance
       rmeas = calc_rmeas(pathname, current_filename, black, white);
       
       % Step 2.2: Correct for substitution error
       rsbse = corr_sbse(rmeas, sbse_coeffs);
       
       % Step 2.3: Calculate corrected erythema index
       eicorr = calc_ei(lambda, rsbse);
       
%        % Step 2.x: Determine uncertainty/weighting?
%        pc_stdev = mua_param(1,:)/5000 + 0.0015; % Percent Standard Deviation
%        Weights = 1./(N.*SE.^2);
%        nonlinmodelW = @(B,t) Weights .* nonlinearmodel(B,t);
%        x = lsqcurvefit(nonlinmodelW,x0,xdata,ydata,lb,ub);
       
       % Step 3: Fit model to measured data
       
       % Step 3.1: Determine coefficients using 'lsqcurvefit'

       [param,~,residual,~,~,~,jacobian] = lsqcurvefit(@(param0,lambda) calc_rd(param0,lambda,mua_param,musp),param0,lambda,rsbse,lb,ub,options);
       %param = real(param);
       %include flag to check if any parameters are at upper or lower bounds?
       
        %rslcf = calc_rd(param,lambda,mua_param,musp);
        %plot(lambda,rsbse,lambda,rslcf)
       
       % Step 3.2: Determine 95% confidence interval using 'nlparci'
       ci = nlparci(param,residual,'jacobian',jacobian);
       paramerror = mean(ci,2) - ci(:,1);
       paramstd = paramerror./1.96;
       
       output(m,:) = [param(1), paramstd(1), param(2), paramstd(2), param(3), paramstd(3), param(4), paramstd(4), param(5), paramstd(5)];
       output(m,11) = eicorr;
       
       m = m+1;
   else
       % Skip this file
   end
   

end

print_to_excel(pathname, xlsfilenames, output)


end

function [pathname, filenames] = select_folder
%SELECT_FOLDER Get a list of files from selected folder

foldername = uigetdir;
pathname = strcat(foldername,'\');
files = dir(foldername);

strucell = struct2cell(files);
strucell1 = strucell(1,:);
filenames = strucell1(3:numel(files));

end

function [lambda, black, white] = get_bgcal(pathname, filenames)
% GET_BGCAL Extracts background and calibration spectra, and wavelength

% Retrieve black spec
fndblk = ones(1,length(filenames)) - double(cellfun('isempty',strfind(filenames, 'black')));
indblk = find(fndblk);
fileblk = strcat(pathname,filenames{indblk});
blkdata = dlmread(fileblk,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
intblk = dlmread(fileblk,' ', [6,3,6,3]); % reads in the integration time, space delimited
blk = (blkdata(:,2)/(intblk/1000))';
black=blk(453:1069);

% Retreive lambda values (500-705)
wl = blkdata(:,1)';
lambda=wl(453:1069);

% Retreive white spec
fndwht = ones(1,length(filenames)) - double(cellfun('isempty',strfind(filenames, 'white')));
indwht = find(fndwht);
filewht = strcat(pathname,filenames{indwht});
whtdata = dlmread(filewht,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
intwht = dlmread(filewht,' ', [6,3,6,3]); % reads in the integration time, space delimited
wht = (whtdata(:,2)/(intwht/1000))';
white=wht(453:1069);

end

function rmeas = calc_rmeas(pathname, current_filename, black, white);
% % CALC_RMEAS Retrieve the modified measured reflectance from Master.Scope
% % files
% 
curfile = strcat(pathname,current_filename);
data = dlmread(curfile,'	', [19,0,2066,1]); % reads in the spectra values, tabs delimited
inttime = dlmread(curfile,' ', [6,3,6,3]); % reads in the integration time, space delimited
spec = (data(:,2)/(inttime/1000))';
spec2 = spec(453:1069);
rmeas = (spec2 - black)./(white - black);

end

function rsbse = corr_sbse(rmeas, sbse_coeffs)

rsbse = ((sbse_coeffs(:,1).*rmeas' + sbse_coeffs(:,2))./(rmeas' + sbse_coeffs(:,3)))';

end

function eicorr = calc_ei(lambda, rsbse)

% Calculate log inverse reflectance
lir = log10(1./rsbse);

% Determine location and values of EI input
lambdaei = [510 543 563 576 610];
for j=1:length(lambdaei)
    [evalue(j), eminloc(j)] = min(abs(lambda-lambdaei(j)));
    eival(j) = mean(lir(eminloc(j)-2:eminloc(j)+2));
end

% Calculate EI
ei = 100*(eival(3) + 1.5*(eival(2) + eival(4)) - 2*(eival(1) + eival(5)));

% Determine location and values of MI input
lambdami = [650 700];
for k=1:length(lambdami)
    [mvalue(k), mminloc(k)] = min(abs(lambda-lambdami(k)));
    mival(k) = mean(lir(mminloc(k)-15:mminloc(k)+15));
end

% Calculate MI
mi = 100*(mival(1)-mival(2) + 0.015);

% Apply MI to get corrected EI
eicorr = ei*(1+0.04*mi);

end

function rslcf = calc_rd(beta,lambda,extco,musp)
%RDSCAT Calculate diffuse reflectance with added scatter losses

% Step 1: Calculate diffuse reflectance using model (rcalc)

% beta(1) = total hemoglobin
% beta(2) = oxygen saturation
% beta(3) = melanin
% beta(4) = background
% beta(5) = background shift

mua = beta(1)*beta(2)*extco(1,:) + beta(1)*(1-beta(2))*extco(2,:) + beta(3)*extco(3,:) + beta(4)*extco(4,:) + beta(5);

mutp = mua + musp;
mueff = sqrt(3.*mua.*mutp);
D = 1./(3.*mutp);

%A = 2.348; % for an nrel = 1.33 (for water)
A = 2.745; % for an nrel = 1.4 (for tissue)

rcalc = musp./((mutp + mueff).*(1 + 2.*A.*D.*mueff));

% Step 2: Convolve spectrum with Gaussian
fwhm = 10;
mu = (lambda(end)-lambda(1))/2 + lambda(1);
goflambda = exp(-(lambda-mu).^2/(2*(fwhm/2.3548)^2));
result = conv(goflambda,rcalc);

cones = ones(1,length(lambda));
cal = conv(goflambda,cones);
norm = result./cal;

rconv = norm(302:918);

% Step 3: Calculate and add in scatter losses

for i=1:length(lambda)
    if mua(i) <= 0.05*musp(i)
        slcf(i) = 0.044*musp(i).^-0.6.*log(mua(i)) + 1.04;
         if slcf(i) >= 1, slcf(i) = 1; end
    else
        slcf(i) = 1;
    end
end

rslcf = slcf.*rconv;

end

function [] = print_to_excel(pathname, xlsfilenames, output)

% Print results to XLSX file
fileparts = strsplit(pathname, '\');
foldername = fileparts(end-1);
xlsname = strcat(pathname, foldername, '.xlsx');
xlsnamestr = xlsname{1};

headings = {'Filename', 'Total Hemoglobin', 'std', 'Oxygen Saturation', 'std', 'Melanin', 'std', 'Background', 'std', 'Shift', 'std', 'EIc'};
xlswrite(xlsnamestr, headings, 'Sheet1', 'A1')

xlsfilenames = xlsfilenames';
xlswrite(xlsnamestr, xlsfilenames, 'Sheet1', 'A2')

xlswrite(xlsnamestr, output, 'Sheet1', 'B2')
end