function [file, pathFile, patchSize, applyTranslationRotation ] = initial_config()
%% initial_config function will initialize all the variables like file name, path & patch size.
%   Keep this file local and not check in as it has machine specific configurations.
applyTranslationRotation = 0;
path = strcat(pwd,'/Multi-Modal-Similarity/Dataset/');
%path = strcat(pwd,'\Dataset\');
file = 'T1_03'; % Change the filename to pick corresponding .TIFF file and it will generate its distance matrix as .dat file
ext = '.TIFF';
fname = strcat(file, ext);
pathFile = strcat(path, fname);
patchSize = 15;

end

