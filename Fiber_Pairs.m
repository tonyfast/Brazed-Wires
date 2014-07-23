%% Fiber Pairs
% Find the fiber pairs

load INDEX/Indexed_fibers.mat

%%
% shift the image to find the neighboring pixels.   there is a function to
% do it, but it is in a toolbox that cant easily be shared.

shifts = [ eye(3); -1 * eye(3) ]

dInoedge = dI( 2:(end-1),2:(end-1),2:(end-1) );

pairs = [0 0];
for ii = 1 : size( shifts, 1 )
    temp = circshift( dI, shifts(ii,:));
    temp = temp( 2:(end-1),2:(end-1),2:(end-1) );
    pairs = union( pairs, unique( [dInoedge(:) temp(:)], 'rows' ) ,'rows');
end

pairs( any(pairs==0,2),:)= []

%%

pairsid = sub2ind( max( dI(:) ) * [1 1], pairs(:,1), pairs(:,2));

G = zeros( max(dI(:)));

G(pairsid) = 1;

%tranpose
pairsid = sub2ind( max( dI(:) ) * [1 1], pairs(:,2), pairs(:,1));

G(pairsid) = 1;

%%

IOI = x0;
reference = dI == IOI;
connected = setdiff( find( G(IOI,:)), IOI );

vol3d( 'Cdata', dI.*ismember( dI,  connected ) + reference )
colormap(([rand(50,3); 0 0 0 ]));

% https://code.google.com/p/yamlmatlab/
units = ReadYaml( 'units.yml' );
daspect( 1./struct2array( units ) )
figure(gcf)
grid on