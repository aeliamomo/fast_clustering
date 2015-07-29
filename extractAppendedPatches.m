function [ patchPairMatrix ] = extractAppendedPatches(isTransformed)
%EXTRACTAPPENDEDPATCHES Summary of this function goes here
%   Detailed explanation goes here
[patchSize,  noOfPosPatches, noOfNegPatches, noOfSample, trainingDataPath, testDataPath] = TuningParametersConfig();

% totalPtsInPatch = patchSize*patchSize;
% pixel_position_x = randi(256);
% pixel_position_y = randi(256);
%noOfSample = noOfNegPatches+1;

base_T1 = strcat(trainingDataPath, 'T1_');
base_T2 = strcat(trainingDataPath, 'T2_');

% base_T1 = 'E:\TUM\Courses\Summer Semester 2015\Machine Learning in Medical Imaging\Project\Extremely Randomized Trees\Dataset\T1_';
% base_T2 = 'E:\TUM\Courses\Summer Semester 2015\Machine Learning in Medical Imaging\Project\Extremely Randomized Trees\Dataset\T2_';

patchPairMatrix = [];
patchPairMatrixT = [];
patchPairMatrixR = [];

translations = [[-60,-60];[-50,-50];[-40,-40];[-30,-30];[-20,-20];[-10,-10];[0,0];[10,10];[20,20];[30,30];[40,40];[50,50];[60,60]];
totalTransformations = size(translations,1);

rotations = [-90;-75;-60;-45;-30;-15;0;15;30;45;60;75;90];


%go through whole 12 images (expect image_8)
for i = 1:12
    %path = strcat(base,num2str(i),'.TIFF');
    if( i ~= 8)
    if(i<10)
        j = strcat('0',num2str(i));
    end
        
    imagePath1 = strcat(base_T1,num2str(j),'.TIFF');
    imagePath2 = strcat(base_T2,num2str(j),'.TIFF');
    
    %get rotated image patch pairs
    boolTranslation = 1;
    [originalImage, translatedImage] = ImageTransformations(imagePath1, translations(i,:), boolTranslation);
    %get rotated image  patche pairs
    boolTranslation = 0;
    [originalImage, rotatedImage] = ImageTransformations(imagePath1, rotations(i), boolTranslation);

if(isTransformed == 0)

for p = 1:noOfPosPatches 
    
    pixel_position_x = randi(256);
    pixel_position_y = randi(256);
    [similarPatches,disSimilarPatches] = extractPatchesPerPixel(imagePath1, imagePath2, pixel_position_x, pixel_position_y, patchSize, noOfSample, 0);
    similarPatches = reshape(cell2mat(similarPatches),[1,2*patchSize*patchSize]); %convert cell to matrix
    disSimilarPatches = reshape(cell2mat(disSimilarPatches),[noOfNegPatches,2*patchSize*patchSize]);%convert cell to matrix
    temp = [similarPatches; disSimilarPatches];
    
    boolAlignedInd = zeros(noOfSample,1);
    boolAlignedInd(1) = 1;
    temp = [temp boolAlignedInd];
    patchPairMatrix = [patchPairMatrix; temp];
    

end
else
    misalignedPatches = extractSimilarPosPatches(originalImage, rotatedImage, patchSize, noOfSample)
    patchPairMatrix = misalignedPatches;
end
    end
end

end

