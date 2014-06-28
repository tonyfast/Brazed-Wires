%%
% Find the center line in the indexed fibers

load INDEX/Indexed_fibers.mat

%% Create Filters
r = 3;
sz = 2*r + 1;


[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;


Rlarge = zeros( size( dI ) );
Rlarge( 1 : sz, 1 : sz, 1 : sz ) = R;
Rlarge(:) = circshift( Rlarge, -1 * ( [r r r] + 0 ) );


R2 =sqrt( rx.^2 + ry.^2 + rz.^2 )<=3;
R2(ceil( numel(R)/2) + 1 ) = 0;


%% Loop over each index

unique_indexes = unique( dI(:) );

for ii = 201%unique_indexes'
    template = simple_convolve(double( dI == ii ), R );
    
    dilated = imdilate( template, R2 );
    
    id = find( abs(template-dilated) >= 1  );
    I = zeros( size(id,1),3);
    [I(:,1), I(:,2), I(:,3)] = ind2sub( size( template ), id );
    ez3dplot( I )
    figure(gcf)
    
    
units = ReadYaml( 'units.yml' );
daspect( 1./struct2array( units ) )

title(numel(id))
    return
end


%%

patterns = [ 1 2 3;
             3 1 2;
             2 3 1 ];

centers = cell( 3, 1);

unid = numel( unique( dI( : ) ) );

for ii = 1 : size( patterns, 1)
    data = permute( dI, patterns( ii, : ) );
    centers{ii} = zeros( unid, size( data,3), 2);
    for zz = 1 : size( data, 3)
        temp = data( :, :, zz );
        [ id ] = find( temp );
        [ ix, iy ] = ind2sub( size( temp ), id );
        centers{ii}(:,zz,1) = accumarray( temp( id ), ix, [unid 1], @mean );
        centers{ii}(:,zz,2) = accumarray( temp( id ), iy, [unid 1], @mean );
        counts{ii}(:,zz) = accumarray( temp( id ), ones(size(ix)), [unid 1], @sum );
    end
end

%% Plot the skeleton of all the fibers


co = cbrewer( 'qual','Set1',3)

for jj = 1 : 3
    id = find(centers{jj}(:,:,1));
    I = zeros( numel( id ), 3 );
    [~, I(:,3)] = ind2sub( size( centers{jj}(:,:,1) ), id );
    temp = centers{jj}(:,:,1);
    I(:,1) = temp( id );
    temp = centers{jj}(:,:,2);
    I(:,2) = temp( id );
    I = permute( I, patterns( jj,: ) );

    for qq = 1 : 3
        b = I(:,qq) > size( dI, qq);
        I(b,:) = [];
    end
    h = ez3dplot( I );
    set( h, 'Markersize', 14, 'MarkerFaceColor',co(jj,:),'MarkerEdgeColor','k');
    if jj == 1;
        hold on;
    end
end
hold off

%%

%% Plot the skeleton of all the fibers


fiberindex = [100: 150];

co = cbrewer( 'qual','Set1',3)
hold on
for jj = 1 : 3
    id = find(centers{jj}(:,:,1));
    
    %%% Filter out fibers
    [ix,~] = ind2sub( size( centers{jj}(:,:,1) ), id );
    
    id = id( ismember( ix, fiberindex ) ); 
    %%%
    
    
    I = zeros( numel( id ), 3 );
    [~, I(:,3)] = ind2sub( size( centers{jj}(:,:,1) ), id );
    temp = centers{jj}(:,:,1);
    I(:,1) = temp( id );
    temp = centers{jj}(:,:,2);
    I(:,2) = temp( id );
    I = permute( I, patterns( jj,: ) );

    for qq = 1 : 3
        b = I(:,qq) > size( dI, qq);
        I(b,:) = [];
    end
    h = ez3dplot( I );
    set( h, 'Markersize', 14, 'MarkerFaceColor',co(jj,:),'MarkerEdgeColor','k');
    if jj == 1;
        hold on;
    end
end
hold off

figure(gcf)