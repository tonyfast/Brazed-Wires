%% Fiber Pairs
% Find the fiber pairs

load INDEX/Indexed_fibers_circular.mat

%%
% shift the image to find the neighboring pixels.   there is a function to
% do it, but it is in a toolbox that cant easily be shared.
%
% IMPROVEMENT: Simple BOOLEAN tests only tell if fibers are touching.
% Change the code to add the number of times the fibers touch.

shifts = [ eye(3); -1 * eye(3) ]

% Remove the edge from the image
dInoedge = dI( 2:(end-1),2:(end-1),2:(end-1) );

pairs = [0 0];

% Fiber index cumulative adjacency matrix
G = zeros( max(dI(:)));
for ii = 1 : size( shifts, 1 )
    disp(sprintf('Shift #%i of #%i.', ii, size(shifts,1)));
    temp = circshift( dI, shifts(ii,:));
    temp = temp( 2:(end-1),2:(end-1),2:(end-1) );
    pairs = [dInoedge(:) temp(:)];
    
    pairs( any(pairs==0,2),:)= [];
    
    G(:) = accumarray( ...
                        sub2ind( ....
                        max( dI(:) ) * [1 1], ....
                        pairs(:,2), pairs(:,1)), ...
                        ones( size(pairs,1), 1), ...
                        [ numel(G), 1], ...
                        @sum );
    
end


% The loop accounds for like-index neighbors.  This statement removes that.
G(:) = G - diag(diag(G));

%% Plot some slices of connected fibers

cut = 300;
[x,y] = find( G>cut, 8, 'first');

FilterIndex = dI .* ismember( dI, ...
    [x(:);y(:)] );


[ trimdex, trimmed] = deal(zeros( size( FilterIndex )));
trimdex = dI;  % Index of the trimmed image
for ii = 1 : numel( x );
    trimmed(:) = max( trimmed ,ii .*[ ismember( FilterIndex, [ x(ii), y(ii) ] )] );
end
% Delete zero slices for a more compact viz
for ii = 1 : 3
    b = all( trimmed == 0,ii);
    for jj = setdiff( 1 : 3, ii )
        deletecol = setdiff( 1 : 3, [ii,jj] );
        b2 = all( b, jj );
        for qq = 1 : 3; evalindex{qq} = ':';end
        evalindex{deletecol} = 'b2';
        evalexp = sprintf( '%s(%s,%s,%s)= [];','trimmed', ...
                            evalindex{1}, evalindex{2}, evalindex{3} );
        eval( evalexp );                
        
        evalexp = sprintf( '%s(%s,%s,%s)= [];','trimdex', ...
            evalindex{1}, evalindex{2}, evalindex{3} );
        eval( evalexp );
        
    end
end
%%

clf;
ax(1) = subplot( 1,2,1)
vol3d( 'Cdata', ...
    trimmed, ...  
    'Alpha', trimmed >0 ...
     );
 
 colormap( rand(10000,3))
 figure(gcf);
 axis equal 
 grid on

 
ax(2) = subplot( 1,2,2)
vol3d( 'Cdata', ...
    trimdex, ...  
    'Alpha', trimmed >0 ...
     );
 
 colormap( rand(10000,3))
 figure(gcf);
 axis equal 
 grid on
 
 linkaxes( ax )

 % dI .* ismember( dI - volumetric image mapped with salient indices
 % [x(:);y(:)] - member set of fibers with # of neighbors > ``cut``
%% Add cumulative adjacency matrix to Indexed fibers

save ./INDEX/Indexed_fibers_circular.mat G -append
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