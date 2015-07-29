function [] = ImageHistogram( input_args )
%IMAGEHISTOGRAM creates and displays histogram of an image based on the
%intensity
I = imread('T2_05.TIFF')
imhist(I)
end

