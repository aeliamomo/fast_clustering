function [I, transformedImage] = ImageTransformations(imagePath, transformation, boolTranslation)
%% Testing - Transformations on an image.
% Translation and Rotation

I=imread(imagePath);
%imshow(I);
  
% Translating an Image
if(boolTranslation)
    transformedImage = imtranslate(I,transformation,'FillValues',0);
    
else % Rotatiting an Image
    transformedImage = imrotate(I,transformation,'crop');
end
%imshow(transformedImage);

end
