% Blend edited and original BedMachine bed topography 

% Load data
load 'DataFiles/bed.mat';
load 'DataFiles/bed_ed.mat';

%Interpolate bed data onto grid
x = min(md.mesh.x):300:max(md.mesh.x);
y = min(md.mesh.y):300:max(md.mesh.y);
bed    = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,bed,x,y,NaN);
bed_ed = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,bed_ed,x,y,NaN);
ocean_mask = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ocean_levelset,x,y,NaN);
ice_mask = InterpFromMeshToGrid(md.mesh.elements,md.mesh.x,md.mesh.y,md.mask.ice_levelset,x,y,NaN);

%Interpolate
distances = zeros(size(x));
transition = 5000;

%Approximate resolution
resolution=x(1,2)-x(1,1);

%set interpolation mask
A = zeros(size(bed));
pos = find(bed~=bed_ed);
A(pos)=1;

%Get distance field
if 1,
   distances1 = resolution*bwdist(A);
   distances2 = resolution*bwdist(~A);
   A= distances1<transition;
   distances1 = resolution*bwdist(~A);
   A= distances2<transition;
   distances2 = resolution*bwdist(~A);
else
   distances1 = resolution*bwdist(~A);
   A= distances1<transition;
   distances2 = resolution*bwdist(~A);
end

%Initialize output
bed_final = bed;

%Linear interpolation in between
pos = find(distances2>0 & distances1>0);
bed_final(pos) = (bed_final(pos).*distances2(pos) + bed_ed(pos).*distances1(pos))./(distances1(pos) + distances2(pos));
bed_final      = single(bed_final);

%Interpolate MITgcm bathymetry exactly
pos = find(distances2==0);
bed_final(pos) = single(bed_ed(pos));

%Save bed and make model consistent
load 'DataFiles/bed.mat';
md.geometry.bed = InterpFromGridToMesh(x',y',bed_final,md.mesh.x,md.mesh.y,0);
pos = find(md.mesh.vertexonboundary);
md.geometry.bed(pos) = bed(pos); clear bed;
pos = find(md.mask.ocean_levelset>0);
md.geometry.base(pos) = md.geometry.bed(pos);
md.geometry.thickness = md.geometry.surface-md.geometry.base;
pos = find(md.geometry.base-md.geometry.bed<0);
md.geometry.bed(pos) = md.geometry.base(pos);

%No negative thickness
pos = find(md.geometry.thickness<=0);
md.geometry.thickness(pos)=1; md.geometry.base(pos) = md.geometry.surface(pos)-md.geometry.thickness(pos);
md.geometry.bed(pos) = md.geometry.base(pos);
pos = find(md.mask.ice_levelset<0 & md.mask.ocean_levelset>0);
md.geometry.base(pos) = md.geometry.bed(pos);
