function [patchSize, noOfPosPatches, noOfNegPatches, noOfSample, trainingDataPath, testDataPath] = TuningParametersConfig()
%TUNING PARAMETERS

patchSize = 15;
%treeDepth = 4;
%noTreeNodes = 2^(treeDepth + 1) - 1;
%totalTreesInForest = 10;

noOfPosPatches = 1;
noOfNegPatches = 10;
noOfSample = noOfNegPatches + 1;

trainingDataPath = strcat(pwd,'/Multi-Modal-Similarity/Dataset/');
testDataPath = strcat(pwd,'/Multi-Modal-Similarity/Dataset/');

end