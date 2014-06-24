%% Braze Wire Data Compaction
% This dataset was shared with by Amanada Levinson at ONR.
% The data is segmented into three phases:
% 
% * Wire
% * Matrix
% * Braze
%
%% Notes on the data
% 
% * Data is versy sparse
% * The braze connects the wires
%
%% Feature identifiers
% 
% # Index and class each fiber.
%% Code objective
% This code converts the tif stack into a volumetric image.

%% File Structure of the Layers 
% Each Layer spacing is 10 microns and each pixel in plane in 5.15 microns

data_dir = 'Braze-Wire-Stack/'; 
to_dir = 'Braze-Wire-Stack-mat/';
unval = 0;
layersfile = dir(fullfile( data_dir, '*.tif'));
ct = 0;

for ii = 1 :1: numel( layersfile )
    ct = ct + 1;
    A = imread( fullfile( data_dir, layersfile(ii).name ) );
    %% reindex the pixel values
    [~, A(:)] = ismember( A, unique(A(:)));
    
    
    if ct == 1
        V = zeros( [ size(A), numel( layersfile)]);
    end
    
    V(:,:,ii) = A;
    
    if mod(ct,20) == 0
        disp( sprintf('Adding layer %i to the stack of %i layers',ii,numel(layersfile)))
    end
    % 3 is the braze
    % 2 is the wire
    
end

return

%% Partition in 2-direction
% Most materials science interactions are short range.  I am going to
% partition the wires to make them easier to save.  The fibers ahve longer
% range on in the 2-direction.  This code block partitions the data for
% easier use.

ct = 0;
for ii = 0 : 500 : size(V,2)
    ct = ct + 1;
    V2 = V(:,(1+ii):min(500+ii,size( V,2)),:);
    save( sprintf('./Braze-Wire-Stack-mat/Braze-Partition-%i-%i.mat',ii+1,ii+500), 'V2')
end