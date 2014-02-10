function [] = model_fit()
%MODEL_FIT Summary of this function goes here
%   Employs fitting algorithm and 1D DRS model to extract chromophore
%   concentrations from total diffuse reflectance spectrum from skin

[pathname, filenames] = select_folder;


end

function [pathcat, strucell2] = select_folder

folder_name = uigetdir;
pathcat = strcat(folder_name,'\');
files = dir(folder_name);

strucell = struct2cell(files);
strucell1 = strucell(1,:);
strucell2 = strucell1(3:numel(files));

end