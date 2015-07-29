% %function [translatedImage,rotatedImage]= translationAndRotation(imagePath,x,y,degree)
% 
% clc;
% clear all;
% close all;
% [patchSize, noOfPosPatches, noOfNegPatches, noOfSample, trainingDataPath, testDataPath] = TuningParametersConfig;
% imagePath = trainingDataPath;
% I = strcat(trainingDataPath,'T1_01.TIFF')
% A=imread(I);
% imshow(I);
% translatedImage = imtranslate(A,[20, 45 ]);
% display('Translation:');
% figure;
% imshow(translatedImage);
% 
% display('Rotation:');
% %using crop to keep image size as original one
% rotatedImage = imrotate(A,135,'crop');
% figure;
% imshow(rotatedImage);
% 
% x = 180;
% y = 50;
% filter = 25;
% halfFilter = floor(filter/2);
% B= A(x-halfFilter:x+halfFilter,y-halfFilter:y+halfFilter);
% size(B)
% imshow(B)
% %rect = getrect


XTest = [];
TestGroundTruthSimilarity = [];

% Perform Transformations on Image for Testing.
translations = [[-60,-60];[-50,-50];[-40,-40];[-30,-30];[-20,-20];[-10,-10];[0,0];[10,10];[20,20];[30,30];[40,40];[50,50];[60,60]];
totalTransformations = size(translations,1);

rotations = [-90;-75;-60;-45;-30;-15;0;15;30;45;60;75;90];

translationSimilarity = [];
rotationSimilarity = [];

for t=1:totalTransformations
    % Translation of Image
    boolTranslationRotation = true;
    [XTransTest, transGroundTruthSimilarity] = SimilarityTestImage(base_T1, patchSize, noPatches, noOfPosPatches, noOfSample, treeDepth, noTreeNodes, totalTreesInForest, structForest, boolTranslationRotation, translations(t,:));
    XTest = [XTest; XTransTest];
    TestGroundTruthSimilarity = [TestGroundTruthSimilarity; transGroundTruthSimilarity];
    
    predTransSimilarity = XTransTest * Weights';
    % Normalise over the number of trees such that the final similarity value is between 0 and 1.
    predNormTransSimilarity = predTransSimilarity./(2^treeDepth);
    meanTransSimilarity = mean(predNormTransSimilarity);
    translationSimilarity = [translationSimilarity; meanTransSimilarity];

    
    % Rotation of Image
    boolTranslationRotation = false;
    [XRotTest, rotGroundTruthSimilarity] = SimilarityTestImage(base_T1, patchSize, noPatches, noOfPosPatches, noOfSample, treeDepth, noTreeNodes, totalTreesInForest, structForest, boolTranslationRotation, rotations(t));
    XTest = [XTest; XRotTest];
    TestGroundTruthSimilarity = [TestGroundTruthSimilarity; rotGroundTruthSimilarity];
    
    predRotSimilarity = XRotTest * Weights';
    % Normalise over the number of trees such that the final similarity value is between 0 and 1.
    predNormRotSimilarity = predRotSimilarity./(2^treeDepth);
    meanRotSimilarity = mean(predNormRotSimilarity);
    rotationSimilarity = [rotationSimilarity; meanRotSimilarity];
end

%% Capture-Range Plot
figure(1);
plot(translations(:,1),translationSimilarity);
%title('Translational-Similarity Capture-Range Plot');
xlabel('Translations');
ylabel('0 \leq Similarity \leq 1');

figure(2);
plot(rotations,rotationSimilarity);
%title('Rotational-Similarity Capture-Range Plot');
xlabel('Rotations');
ylabel('0 \leq Similarity \leq 1');

predictedSimilarity = XTest * Weights';

% Normalise over the number of trees such that the final similarity value is between 0 and 1.
predNormSimilarity = predictedSimilarity./(2^treeDepth);


%% Evaluation - Performance Measures!

%% =========== Classification - Confusion Matrix =============
% Define the vectors for the ground truth and the predicted class for each sample.
groundTruth = TestGroundTruthSimilarity; %[1,0,0,1,0,0,1,1,0,1,0,0,1,1,0,0,1,1,0,1,0,1]';
Predictions = predNormSimilarity; %predictedSimilarity; %[0.6,0.2,0.3,0.5,0.9,0.8,0.6,0.3,0.1,0.2,0.1,0.5,0.8,0.1,0.3,0.7,0.8,0.2,0.4,0.4,0.5,0.4]';

thresholdVec = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
totalThresholds = size(thresholdVec,2);

lstconfMat2by2 = {};
sensitivityVec = zeros(1, totalThresholds);
specificityVec = zeros(1, totalThresholds);
PPVVec = zeros(1, totalThresholds);
NPVVec = zeros(1, totalThresholds);
accuracyVec = zeros(1, totalThresholds);
F1_MeasureVec = zeros(1, totalThresholds);

for t=1:totalThresholds
    
    % =========== Part a: Threshold Confusion Matrix =============
    % The function threshold_confusion_matrix expects a vector containing the ground truth and
    % a vector of containing the predicted probabilities for each sample.
    % For each unique threshold the function should return a 2 ? 2 confusion matrix.

    confMat2by2 = threshold_confusion_matrix(groundTruth, Predictions, thresholdVec(t));
    lstconfMat2by2 = [lstconfMat2by2 ; confMat2by2];
    
    % =========== Part b: Performance Measures =============
    % For any number of classes calculate the following 6 measures based on a confusion matrix:
    % Sensitivity, specificity, positive predictive value, negative predictive value, accuracy, and F1 measure.
    [sensitivityVec(t), specificityVec(t), PPVVec(t), NPVVec(t), accuracyVec(t), F1_MeasureVec(t)] = performanceMeasures(lstconfMat2by2{t});
end


%% =========== ROC and Precision-Recall Curve =============
% For a binary classifier - Let the classes be represented by 0 and 1.

% Part a. ROC
% False positive rate = 1 - Specificity
figure(3);
plot(1 - specificityVec, sensitivityVec);
title('ROC');

% Part b. Precision-Recall curve
% Precision = Positive predictive value (PPV)
% Recall = Sensitivity !!
figure(4);
plot(PPVVec, sensitivityVec);
title('Precision-Recall curve');






%end