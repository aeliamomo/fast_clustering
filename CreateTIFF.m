function [] = CreateTIFF()
%Reads all values from the directory and creates TIFF images. We have the
%ability to set a few HEADERS in the TIFF images we generate but is is
%currently disabled. The setting of he HEADERS would give us the ability to
%control the number of STRIPS that our Block Extraction Algorithm reads at
%a time. USAGE: CreateTIFF

% ===== NOTE: PLease modify the directory of read images below, as per needs ======
d=dir(fullfile(pwd,'/Multi-Modal-Similarity/Dataset/*.nii'));
for k=1:length(d)
  fname=d(k).name;
  path = strcat(pwd,'/Multi-Modal-Similarity/Dataset/');
  pathFile = strcat(path,fname);
  Data = load_untouch_nii(pathFile);
  %imshow(double(Data.img)/255);
  WriteFname = strcat(path,fname(1:end-4),'.TIFF');
  
  %%set TIFF specific HEADERS-- DISABLED for now.
  %tiffHandle = Tiff(WriteFname,'w');
  %tagstruct.ImageLength = size(single(Data.img)/255,1)
  %tagstruct.ImageWidth = size(single(Data.img)/255,2)
  %tagstruct.Photometric = Tiff.Photometric.MinIsBlack
  %tagstruct.BitsPerSample = 32
  %tagstruct.SamplesPerPixel = 1
  %tagstruct.RowsPerStrip = 2
  %tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
  %tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky
  %tagstruct.Software = 'MATLAB'
  %tiffHandle.setTag(tagstruct)
  %tiffHandle.write(single(Data.img)/255);
  
  %%write the TIFF images 
  imwrite(double(Data.img)/255,WriteFname);
end
    
end