function point_cloud = gen_pointcloud(mat3d,thresh)
if ~isreal(mat3d)
    mat3d = abs(mat3d);
end
if nargin<2
    thresh = ( max(mat3d,[],'all')+min(mat3d,[],'all') )/4;
end

% 3D input:
[row,col,vol] = ind2sub(size(mat3d),find(mat3d > thresh));
pixel = zeros(length(row),1);
for ii = 1:length(row)
    pixel(ii,1)=mat3d(row(ii),col(ii),vol(ii));
end
point_cloud = [row,col,vol,pixel];


% 4D input:
[row,col,vol1,vol2] = ind2sub(size(mat3d),find(mat3d > thresh));
pixel = zeros(length(row),1);
for ii = 1:length(row)
    pixel(ii,1)=mat3d(row(ii),col(ii),vol1(ii),vol2(ii));
end
point_cloud = [row,col,vol1,vol2,pixel];

% 6D input:
[d1,d2,d3,d4,d5,d6] = ind2sub(size(mat3d),find(mat3d > thresh));
pixel = zeros(length(d1),1);
for ii = 1:length(d1)
    pixel(ii,1)=mat3d(d1(ii),d2(ii),d3(ii),d4(ii),d5(ii),d6(ii));
end
point_cloud = [d1,d2,d3,d4,d5,d6,pixel];

end