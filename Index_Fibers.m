if ~exist('V2','var')
%     load Braze-Wire-Stack-mat/Braze-Partition-1501-2000.mat
    load Braze-Wire-Stack-mat/CompleteStack.mat
end

%%  Binarize Image

inc = 3;

A = V2(1:inc:end,1:inc:end,1:inc:end) == 2;

% %% Spherical Filter
% 
% % filter radius
% r = 5;
% sz = 2*r + 1;
% 
% [rx, ry, rz] = meshgrid( -r:r );
% R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;
% 
% Rlarge = zeros( size( A ) );
% Rlarge( 1 : sz, 1 : sz, 1 : sz ) = R;
% Rlarge(:) = circshift( Rlarge, -1 * ( [r r r] + 0 ) );
% 
% %% Use SpatialStatistics as Spherical Hough Filter
% tic;
% disp('Started spatial statistics..')
% Filter = SpatialStatsFFT( A, Rlarge, 'display', false, ...
%                                      'normalize', false);
% toc;
% disp( sprintf('Filter max : %f', max(Filter(:))));
% disp( sprintf('Area of sphere: %f', sum(R(:) )));
% 
% %%  Find centers of fibers using 2:D filters
% 
% Rlarge(:) = 0;
% [rx, ry] = meshgrid( -r:r );
% R = sqrt( rx.^2 + ry.^2  )<=r;
% Rlarge( 1 : sz, 1 : sz ) = R;
% Rlarge(:) = circshift( Rlarge, -1 * ( [r r ] + 0 ) );
% Rlarge(1) = 0;
% 
% patterns = [ 1 2 3;
%              1 3 2;
%              3 1 2 ];
% 
% I = cell( size( patterns,2 ),1);
% for pp = 1 : size( patterns,2 )
%     % Find peaks
%     disp('Finding Peaks')
%     tic;
%     Stats = imdilate( Filter, Rlarge);
%     toc;
%     
%     id = find( Stats < Filter );
%     I{pp} = zeros( size(id,1),3);
%     
%     [I{pp}(:,1), I{pp}(:,2), I{pp}(:,3)] = ind2sub( size( Filter ), id );
% end
%              
% 
% 
% %% 3-D Filter
% R2 = R;
% R2(ceil( numel(R)/2) + 1 ) = 0;
% dilated = imdilate( Filter, R2 );
% 
% id = find( Filter >= dilated );
% I = zeros( size(id,1),3);
% [I(:,1), I(:,2), I(:,3)] = ind2sub( size( Filter ), id );
% 

%%

r = 2.7;
sz = 2*r + 1;

[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;


fA = round(simple_convolve( A, R ));

thres_fA = fA.*A > .9*sum(R(:));

clf;
vol3d('Cdata',thres_fA)
figure(gcf);

%% Eroded Fibers
% Erode the original fibers to index the middle then watershed the values
% out.

r = 2;
sz = 2*r + 1;

[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<=r;

% Initialize with filtered image
A2 = A;
A2 = thres_fA;
eroded = imerode( A2, R);
id = find( eroded );
I = zeros( size(id,1),3);
[I(:,1), I(:,2), I(:,3)] = ind2sub( size( A2 ), id );

L = bwlabeln( eroded );

clf; 
vol3d('Cdata',L);
colormap( rand(1000,3));
figure(gcf)

%% Remove low populated indexes

index_count = accumarray( L(L~=0), ones( sum(L(:)~=0),1),[],@sum)

ct_cut = 100;
low_id = find( index_count < ct_cut );

L2 = L;
% reindex labelled image
L2(ismember(L2,low_id)) = 0;


%% Plot eroded data
% clf;
% hold on
% for ii = unique( L(L(:)>0) )'
%     id = find( L == ii );
%     I = zeros( size(id,1),3);
%     [I(:,1), I(:,2), I(:,3)] = ind2sub( size( A2 ), id );
%     
%     h = ez3dplot(I);
%     
%     set( h, 'Color',rand(1,3));
% 
% end
% hold off
% axis equal; 
% grid on
% figure(gcf)

%%  Dilate Fibers to roughly approximate the complete fiber index 
% Incrementally dilate the indices
%
% Rules
% # Do not replace one index with another index.  We'll never know really.

r = 2;
[rx, ry, rz] = meshgrid( -r:r );
R = sqrt( rx.^2 + ry.^2 + rz.^2 )<r;

dI = L2;% dilate fiber indexed image'
f = {@mean, @std, @skewness, @kurtosis}
for ii = 1 : 20
    temp = imdilate( dI, R ) .* A2;  % Add a mask of the original image
    % placeholder for what is happening on initializing
    %     b = temp == dI;
    %     dI(b) = temp(b);
    %     
    % Real work
    % When a new index has been made
    b2 = temp(:) > 0 & dI(:) > 0 & temp(:)~=dI(:);
    id = find(b2);
    [x y z] = ind2sub(size(dI), id );
%     
%     ida = find(dI(:));
%     [xa ya za] = ind2sub(size(dI), ida );
%     % Choose new label based on minimal contribute to the std
%     % Baseline for prior image
%     for ff = 1 : numel( f)
%         mx(:,ff) = accumarray( dI(id), x, [ max(L(:)) 1], f{ff});
%         my(:,ff) = accumarray( dI(id), y, [ max(L(:)) 1], f{ff});
%         mz(:,ff) = accumarray( dI(id), z, [ max(L(:)) 1], f{ff});
%     end
%     
%     for ff = 1 : numel( f)
%         mxa(:,ff) = accumarray( dI(ida), xa, [ max(L(:)) 1], f{ff});
%         mya(:,ff) = accumarray( dI(ida), ya, [ max(L(:)) 1], f{ff});
%         mza(:,ff) = accumarray( dI(ida), za, [ max(L(:)) 1], f{ff});
%     end
%     
%     [mx(isnan(mx)), my(isnan(my)), mz(isnan(mz))] = deal(0);
    
%     if any(b2)  
%         break;
%     end
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
colormap( rand(1000,3));
figure(gcf);
axis equal
grid on;

%%
units = ReadYaml( 'units.yml' );
daspect( 1./struct2array( units ) )
%% Units

b = all( mx==0,2) & all(my==0,2) & all(mz==0,2);
