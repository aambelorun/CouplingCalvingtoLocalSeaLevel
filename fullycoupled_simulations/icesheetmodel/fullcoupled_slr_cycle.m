if ~exist('steps', 'var') || isempty(steps)
	steps = [2]; 
end

if ~exist('loadonly', 'var') || isempty(loadonly)
	loadonly = 1;
end

if ~exist('sigma', 'var') || isempty(sigma)
        sigma = 230;
end

if ~exist('k', 'var')
        error('Cycle index k not provided from bash.');
end

codepath = '';%path to directory with ISSM binaries
execpath = '';%path to directory for simulation execution
coresrequested       = 24;
memrequested         = 5; %GB
timerequested        = 10*60; %minutes
%mycluster = pace('np',coresrequested,'login','loginname','codepath',codepath,'executionpath',execpath,'mem',memrequested,'time',timerequested);
mycluster = generic('name', oshostname(),'np',4);

org=organizer('repository','ModelsSeaTest','prefix', sprintf('fullSlr_%d_', sigma), 'steps',steps);
addpath('../../Functions');

% {{{
if perform(org,sprintf('fromdh%d', k))

    mds = loadmodel(sprintf('../sealevelmodel/Models/fullSlr_%d_THICKNESS%d_Transient', sigma, k));
    md = loadmodel(org, sprintf('fromdh%d', k-1));

    %%% Set initial conditions %%%
    md.geometry.thickness        = md.results.TransientSolution(end).Thickness;
    md.initialization.vx         = md.results.TransientSolution(end).Vx;
    md.initialization.vy         = md.results.TransientSolution(end).Vy;
    md.initialization.vel        = md.results.TransientSolution(end).Vel;
    md.mask.ocean_levelset       = md.results.TransientSolution(end).MaskOceanLevelset;
    md.initialization.pressure   = md.results.TransientSolution(end).Pressure;
    md.geometry.base             = md.results.TransientSolution(end).Base;
    md.geometry.surface          = md.results.TransientSolution(end).Surface;
    md.mask.ice_levelset         = md.results.TransientSolution(end).MaskIceLevelset;

    sealevel = mds.results.TransientSolution(end).Sealevel-mds.results.TransientSolution(end).Bed;

    disp('Projecting slr to mesh...');

    pos = find(mds.mesh.lat < 0);
    lat_antarc = mds.mesh.lat(pos);
    lon_antarc = mds.mesh.long(pos);
    data_antarc = sealevel(pos);

    disp('Projecting global coordinates...');
    [x_south, y_south] = ll2xy(lat_antarc, lon_antarc, -1);
    elems_orig = mds.mesh.elements;
    mask = all(ismember(elems_orig, pos), 2);
    elems_kept = elems_orig(mask, :);

    [~, elems_new] = ismember(elems_kept, pos);

    disp('Interpolating sea level');
    md.initialization.sealevel = InterpFromMeshToMesh2d(elems_new, x_south, y_south, data_antarc, md.mesh.x, md.mesh.y, 'default', 0);

    md.calving=calvingvonmises();
    md.calving.stress_threshold_groundedice = 1e6;
    md.calving.stress_threshold_floatingice = sigma*1e3;

    md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
    md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
    md.levelset.migration_max = 1000000.0;
    md.levelset.stabilization = 1;
    md.transient.ismovingfront = 1;
    md.transient.isgroundingline = 1;
    md.timestepping.time_step = 1/12;
    md.timestepping.final_time = 1;
    md.settings.output_frequency=1;
    md.levelset.reinit_frequency = 20;

    md.miscellaneous.name = sprintf('fullSlr_%d_fromdh%d', sigma, k);
    md.settings.waitonlock = 0;
    md.verbose=verbose('solution',true,'convergence',true);
    md.cluster = mycluster;
    md.transient.requested_outputs={'default','IceVolume','GroundinglineMassFlux','IceVolumeAboveFloatation','TotalCalvingFluxLevelset'}; %For rate-based calving laws
    md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

    if loadonly
		%md=loadresultsfromcluster(md);
        savemodel(org,md);
    end
end
%}}}

