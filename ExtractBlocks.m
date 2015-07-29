function patches = ExtractBlocks(filename, patchSize)
%EXTRACTBLOCKS extracts the blocks from the images and returns patches in a
% matrix. whose each row corresponds to a patch.
%USAGE: ExtractBlocks('T1_01.TIFF')

%%extract blocks and reshape them to flat out array
%Changes to make patchSize dynamic by passing it as an input parameter.
%extract = @(block_struct) reshape(block_struct.data,[1 100]);

totalPointsInPatch = patchSize * patchSize;
extract = @(block_struct) reshape(block_struct.data,[1 totalPointsInPatch]);

%Changes to make patchSize dynamic by passing it as an input parameter.
%Blocks = blockproc(filename,[10 10],extract,'PadPartialBlocks',true);

Blocks = blockproc(filename,[patchSize patchSize],extract,'PadPartialBlocks',true);

%   for i= 1:16 
%     imshow(reshape(Blocks(i,256*8+1:256*9),[16,16]));
%     hold on;
%   end

%%blocproc concats into a matrix, so flat out the matrix

pixels = numel(Blocks);

%Blocks = reshape(Blocks,[1,pixels])

[row,col] = size(Blocks);

X = [];

for i = 1:row

  for j = 1:row

      X = [X ,Blocks(i,(j-1)*totalPointsInPatch+1 :j*totalPointsInPatch)];

  end

end

%   for i = 1:row
%       X = [X ,Blocks(1,(i-1)*totalPointsInPatch+1 :i*totalPointsInPatch)];
%   end 
%Blocks = reshape(Blocks,[1 pixels]);
%%extract and arrange the patch specific pixels in a row
%   patches = zeros(pixels/100,100);
%   offset=1;
%   for i=1:(pixels/100)
%       patches(i,:)=Blocks(offset:(offset+99));
%       offset =offset+100;
%   end

Blocks = X;
pixels = size(Blocks,2);
%Changes to make patchSize dynamic by passing it as an input parameter.

%%extract and arrange the patch specific pixels in a row

patches = zeros(pixels/totalPointsInPatch,totalPointsInPatch);
offset=1;
for i=1:(pixels/totalPointsInPatch)
  patches(i,:)=Blocks(offset:(offset+totalPointsInPatch-1));
  offset =offset+totalPointsInPatch;
end

end