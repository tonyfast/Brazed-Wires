if ~exist('V2','var')
%     load Braze-Wire-Stack-mat/Braze-Partition-1501-2000.mat
    load Braze-Wire-Stack-mat/CompleteStack.mat
end

%%  Binarize Image

inc = 3;

A = V2(1:inc:end,1:inc:end,1:inc:end) == 2;

%% Spherical Filter

% filter radius
r = 5;
sz = 2*r + 1;

[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;

Rlarge = zeros( size( A ) );
Rlarge( 1 : sz, 1 : sz, 1 : sz ) = R;
Rlarge(:) = circshift( Rlarge, -1 * ( [r r r] + 0 ) );

%% Use SpatialStatistics as Spherical Hough Filter
tic;
disp('Started spatial statistics..')
Filter = SpatialStatsFFT( A, Rlarge, 'display', false, ...
                                     'normalize', false);
toc;
disp( sprintf('Filter max : %f', max(Filter(:))));
disp( sprintf('Area of sphere: %f', sum(R(:) )));

%%  Find centers of fibers using 2:D filters

Rlarge(:) = 0;
[rx, ry] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2  )<=r;
Rlarge( 1 : sz, 1 : sz ) = R;
Rlarge(:) = circshift( Rlarge, -1 * ( [r r ] + 0 ) );
Rlarge(1) = 0;

patterns = [ 1 2 3;
             1 3 2;
             3 1 2 ];

I = cell( size( patterns,2 ),1);
for pp = 1 : size( patterns,2 )
    % Find peaks
    disp('Finding Peaks')
    tic;
    Stats = imdilate( Filter, Rlarge);
    toc;
    
    id = find( Stats < Filter );
    I{pp} = zeros( size(id,1),3);
    
    [I{pp}(:,1), I{pp}(:,2), I{pp}(:,3)] = ind2sub( size( Filter ), id );
end
             


%% 3-D Filter
R2 = R;
R2(ceil( numel(R)/2) + 1 ) = 0;
dilated = imdilate( Filter, R2 );

id = find( Filter >= dilated );
I = zeros( size(id,1),3);
[I(:,1), I(:,2), I(:,3)] = ind2sub( size( Filter ), id );

%% Eroded Fibers
% Erode the original fibers to index the middle then watershed the values
% out.

A2 = A;
eroded = imerode( A2, ones(9));
id = find( eroded );
I = zeros( size(id,1),3);
[I(:,1), I(:,2), I(:,3)] = ind2sub( size( A2 ), id );

L = bwlabeln( eroded );

clf;
hold on
for ii = unique( L(L(:)>0) )'
    id = find( L == ii );
    I = zeros( size(id,1),3);
    [I(:,1), I(:,2), I(:,3)] = ind2sub( size( A2 ), id );
    
    h = ez3dplot(I)
    
    set( h, 'Color',rand(1,3));

end
hold off
axis equal; 
grid on
figure(gcf)

%%  Dilate Fibers to roughly approximate the complete fiber index 
% Incrementally dilate the indices
%
% Rules
% # Do not replace one index with another index.  We'll never know really.

r = 2;
[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;

dI = L;% dilate fiber indexed image
for ii = 1 : 20
    temp = imdilate( dI, R ) .* A2;  % Add a mask of the original image
    % placeholder for what is happening on initializing
    %     b = temp == dI;
    %     dI(b) = temp(b);
    %     
    % Real work
    % When a new index has been made
    b = temp > 0 & dI == 0;
    dI(b) = temp(b);
    %     % Conflicting indices
    %     b = temp >0 & temp ~= dI;
    %     dI(b) = dI(b)
    
       % Assign the positions that are already good
    % overlap of indices
    disp( ii );
end

%%

clf;
vol3d( 'Cdata', dI );
colormap( rand(100,3));
figure(gcf);
axis equal
grid on;


units = ReadYaml( 'units.yml' );
daspect( 1./struct2array( units ) )
%% Units

