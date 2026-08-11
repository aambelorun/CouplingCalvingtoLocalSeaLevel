if ~exist('steps', 'var') || isempty(steps)
	steps = [1]; 
end

if ~exist('loadonly', 'var') || isempty(loadonly)
	loadonly = 1;
end

if ~exist('k', 'var') || isempty(k)
        k=1;
end

if ~exist('cf', 'var') || isempty(cf)
        cf = 5;
end

if ~exist('gl', 'var') || isempty(gl)
        gl = 0;
end

if ~exist('sigma', 'var') || isempty(sigma)
        sigma = 230;
end

codepath = '';%path to directory with ISSM binaries
execpath = '';%path to directory for simulation execution
coresrequested       = 24;
memrequested         = 5; %GB
timerequested        = 10*60; %minutes
%mycluster = pace('np',coresrequested,'login','loginname','codepath',codepath,'executionpath',execpath,'mem',memrequested,'time',timerequested);
mycluster = generic('name', oshostname(),'np',4);

expname = sprintf('cf%dgl%d', cf, gl);

org=organizer('repository','Models','prefix', sprintf('idealized_%d_', sigma), 'steps',steps);
addpath('../Functions');

% {{{

if perform(org, sprintf('%s_initial', expname))
org
	md = loadmodel('../inputs/initialization/Models/Amundsen_Transient.mat');

   	if cf > 0 && gl == 0

        flags = sealevel_mask(md,1e3);
        md.initialization.sealevel = zeros(md.mesh.numberofvertices,1);
        md.initialization.sealevel(flags) = -cf;

    elseif cf == 0 && gl > 0
        flags = sealevel_mask(md,1e3);
        md.initialization.sealevel = -gl * ones(md.mesh.numberofvertices,1);
        md.initialization.sealevel(flags) = 0;
   	end 
   
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

   	md.miscellaneous.name = sprintf('idealized_%d_%s_1', sigma, expname);
   	md.settings.waitonlock = 0;
   	md.verbose=verbose('solution',true,'convergence',true);
   	md.cluster = mycluster;
   	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'}; 
   	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

   	if loadonly
    	%md=loadresultsfromcluster(md);
    	savemodel(org,md);
   	end
end
%}}}

% {{{
if perform(org,sprintf('%s_%d',expname, k))

    if k==2
		md = loadmodel(org, sprintf('%s_initial', expname));
	else
		md = loadmodel(org, sprintf('%s_%d', expname, k-1));	
	end
	
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
        
	if cf > 0 && gl == 0
		flags = sealevel_mask(md,1e3);
        md.initialization.sealevel = zeros(md.mesh.numberofvertices,1);
		md.initialization.sealevel(flags) = -cf;

	elseif cf == 0 && gl > 0
        flags = sealevel_mask(md,1e3);
        md.initialization.sealevel = -gl * ones(md.mesh.numberofvertices,1);
        md.initialization.sealevel(flags) = 0;
	end

	    md.calving=calvingvonmises();
        md.calving.stress_threshold_groundedice = 1e6;
        md.calving.stress_threshold_floatingice = sigma*1e3;

        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.levelset.spclevelset = NaN(md.mesh.numberofvertices,1);
        md.levelset.migration_max = 1000000.0; %avoid too fast advance/retreat of the front
        md.levelset.stabilization = 1;
        md.transient.ismovingfront = 1;
        md.transient.isgroundingline = 1;
        md.timestepping.time_step = 1/12;
        md.timestepping.final_time = 1;
	    md.settings.output_frequency=1;
        md.levelset.reinit_frequency = 20;

        md.miscellaneous.name = sprintf('idealized_%d_%s_%d', sigma, expname, k);
        md.settings.waitonlock = 0;
        md.verbose=verbose('solution',true,'convergence',true);
        md.cluster = mycluster;
        md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'}; 
        md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

        if loadonly
 		    %md=loadresultsfromcluster(md);
            savemodel(org,md);
        end
end
%}}}

% {{{

if perform(org, sprintf('%s', expname))

        md = loadmodel('../inputs/initialization/Models/Amundsen_Transient.mat');

        md.initialization.sealevel = -cf * ones(md.mesh.numberofvertices,1);
	
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
        md.timestepping.final_time = 2;
        md.settings.output_frequency=1;
	    md.levelset.reinit_frequency = 20;

        md.miscellaneous.name = sprintf('idealized_%d_%s', sigma, expname);
        md.settings.waitonlock = 0;
        md.verbose=verbose('solution',true,'convergence',true);
        md.cluster = mycluster;
        md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation'};
        md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

        if loadonly
                %md=loadresultsfromcluster(md);
                savemodel(org,md);
        end
end
%}}}
