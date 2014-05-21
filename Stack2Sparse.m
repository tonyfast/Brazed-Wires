%% Braze Wire Data Compaction
% This dataset was shared with by Amanada Levinson at ONR.
% The data is segmented into three phases:
% 
% * Wire
% * Matrix
% * Braze
%
% The data is Three Phases and very sparse.  I am converting this data into
% smaller sparse matrices.


%% File Structure of the Layers 
% Each Layer spacing is 10 microns and each pixel in plane in 5.15 microns

data_dir = 'Braze-Wire-Stack/'; 
to_dir = 'Braze-Wire-Stack-mat/';
unval = 0;
layersfile = dir(fullfile( data_dir, '*.tif'));
ct = 0;
for ii = 1 :5: numel( layersfile )
    ct = ct + 1;
    A = imread( fullfile( data_dir, layersfile(ii).name ) );
    %% reindex the pixel values
    [~, A(:)] = ismember( A, unique(A(:)));
    % 3 is the braze
    % 2 is the wire
    
    layer.braze = find( A(:) == 3); 
    layer.wire = find( A(:) == 2);
    layer.size = size(A);
    
    save( fullfile( to_dir, regexprep( layersfile(ii).name, '.tif','.mat')), 'layer');
   return
end