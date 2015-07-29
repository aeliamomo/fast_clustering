function [ output_args ] = patchWiseDiff( input_args )
%PATCHWISEDIFF Summary of this function goes here
%   Detailed explanation goes here

 extract = @(block_struct) reshape(block_struct.data,[1 totalPointsInPatch]);
 Blocks = blockproc('T1_01.TIFF',[100 100],patchRegistration,'PadPartialBlocks',true);
end

