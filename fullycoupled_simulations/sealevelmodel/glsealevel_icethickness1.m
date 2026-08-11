%clear all;
if ~exist('steps', 'var') || isempty(steps)
        steps = [1:3 5]; 
end

if ~exist('loadonly', 'var') || isempty(loadonly)
        loadonly = 0;
end

if ~exist('sigma', 'var') || isempty(sigma)
        sigma =230;
end

codepath = '';%path to directory with ISSM binaries
execpath = '';%path to directory for simulation execution
coresrequested       = 24;
memrequested         = 5; %GB
timerequested        = 10*60; %minutes
%mycluster = pace('np',coresrequested,'login','loginname','codepath',codepath,'executionpath',execpath,'mem',memrequested,'time',timerequested);
mycluster = generic('name', oshostname(),'np',4);

org=organizer('repository','Models','prefix', sprintf('glSlr_%d_THICKNESS1_', sigma),'steps',steps);

% {{{ Mesh:
if perform(org,'Mesh'),
	disp('   Step 1: Global mesh creation');

	numrefine = 1;
	resolution = 150*1e3;		% inital resolution [m]
	radius = 6.371012*10^6;		% mean radius of Earth, m
	mindistance_coast = 25*1e3;	% coastal resolution [m]
	mindistance_land = 300*1e3;	% resolution on the continents [m]
	maxdistance = 600*1e3;		% max element size (on mid-oceans) [m]

	% get resolution from amundsen simulation

	mds = loadmodel('Models/AmundsenSeaTestMesh');
	areas = GetAreas(mds.mesh.elements,mds.mesh.x,mds.mesh.y);
	res_am_elements = sqrt(2.*areas);
	res_am_nodes = averaging(mds, res_am_elements, 0);

	am_center_x = mean(mds.mesh.x);
	am_center_y = mean(mds.mesh.y);
	am_halo_radius = 1000 * 1e3; 
	am_halo_res    = 30 * 1e3;  
    
   	% initial global mesh 
	md=model;
	md.mesh=gmshplanet('radius',radius*1e-3,'resolution',resolution*1e-3); 
   
	for i=1:numrefine, % refine mesh twice

        %coastal refinement
		ocean_mask = gmtmask(md.mesh.lat, md.mesh.long); 
		distance=zeros(md.mesh.numberofvertices,1);

		pos=find(~ocean_mask); % land
    	coaste.lat=md.mesh.lat(pos);
   		coaste.long=md.mesh.long(pos);
   		pos=find(ocean_mask);
   		coasto.lat=md.mesh.lat(pos);
   		coasto.long=md.mesh.long(pos);

    	% coastal distance calculation
		for j=1:md.mesh.numberofvertices
			phi1 = md.mesh.lat(j)/180*pi; lambda1=md.mesh.long(j)/180*pi;
			if ocean_mask(j),
				phi2 = coaste.lat/180*pi; lambda2=coaste.long/180*pi;
				deltaphi=abs(phi2-phi1); deltalambda=abs(lambda2-lambda1);
				d=radius*2*asin(sqrt(sin(deltaphi/2).^2+cos(phi1).*cos(phi2).*sin(deltalambda/2).^2));
			else
				phi2 = coasto.lat/180*pi; lambda2=coasto.long/180*pi;
				deltaphi=abs(phi2-phi1); deltalambda=abs(lambda2-lambda1);
				d=radius*2*asin(sqrt(sin(deltaphi/2).^2+cos(phi1).*cos(phi2).*sin(deltalambda/2).^2));
			end
			distance(j)=min(d);
        end

        distance(distance < mindistance_coast) = mindistance_coast;
        distance(ocean_mask~=1 & distance > mindistance_land) = mindistance_land;
        dist = min(maxdistance, distance);
        
        refine_metric = dist/2;

        [x_global, y_global]= ll2xy(md.mesh.lat, md.mesh.long, -1);
		
        dist_to_am_center = sqrt((x_global - am_center_x).^2 + (y_global - am_center_y).^2);
        in_halo = find(dist_to_am_center < am_halo_radius & md.mesh.lat < 0);
        refine_metric(in_halo) = min(refine_metric(in_halo), am_halo_res);

        % interpolate target resolution from local mesh to global points
        am_res_on_global = InterpFromMeshToMesh2d(mds.mesh.elements, mds.mesh.x,...
        mds.mesh.y, res_am_nodes, x_global, y_global, 'default', NaN);

        inamundsen = find(~isnan(am_res_on_global) & md.mesh.lat < 0);
        refine_metric(inamundsen) = am_res_on_global(inamundsen);

        md.mesh=gmshplanet('radius',radius*1e-3,'resolution',resolution*1e-3,'refine',md.mesh,'refinemetric',refine_metric);
	end

    ocean_mask=gmtmask(md.mesh.lat,md.mesh.long);
	pos = find(ocean_mask==0);
   	md.mask.ocean_levelset=-ones(md.mesh.numberofvertices,1);
   	md.mask.ocean_levelset(pos)=1;
	savemodel(org,md);

	%plotmodel (md,'data',md.mask.ocean_levelset,'edgecolor','k','view',[45 45]);

end % }}}

% {{{Loads:
if perform(org,'Loads'),  
	disp('   Step 2: Define loads in meters of ice height equivalent');
	md = loadmodel(org, 'Mesh');
	mds = loadmodel(sprintf('../icesheetmodel/ModelsSeaTest/glSlr_%d_fromobs', sigma)); %initial 1 year simulation

	[x_global, y_global]= ll2xy(md.mesh.lat, md.mesh.long, -1);

	h0 = mds.geometry.thickness; % Initial thickness
	hend = mds.results.TransientSolution(end).Thickness;
	dh_amundsen = hend - h0;

	ice_load = InterpFromMeshToMesh2d(mds.mesh.elements, mds.mesh.x,...
	mds.mesh.y, dh_amundsen, x_global, y_global, 'default', 0);
	pos_north = find(md.mesh.lat > 0);
	ice_load(pos_north) = 0;
	
	md.geometry.bed=zeros(md.mesh.numberofvertices,1);
	md.geometry.base=md.geometry.bed;
	md.geometry.thickness = 100*ones(md.mesh.numberofvertices,1);
	md.geometry.surface=md.geometry.bed+md.geometry.thickness;

	md.masstransport.spcthickness = repmat(md.geometry.thickness,1,2); 
	md.masstransport.spcthickness(:,end) = md.masstransport.spcthickness(:,end) + ice_load; 
	md.masstransport.spcthickness(end+1,:) = [0 1]; 

	md.smb.mass_balance=zeros(md.mesh.numberofvertices,1); 
	savemodel(org,md);
end % }}}

% {{{ Parameterization:
if perform(org,'Parameterization'),  
	disp('   Step 3: Parameterization');
	md = loadmodel(org,'Loads');
	
	md.mask.ice_levelset=-md.mask.ocean_levelset;
	
	md.solidearth.lovenumbers = lovenumbers('maxdeg',10000);
	md.solidearth.settings.reltol=NaN;
	md.solidearth.settings.abstol=1e-3;
	md.solidearth.settings.sealevelloading=1;
	md.solidearth.settings.isgrd=1;
	md.solidearth.settings.grdmodel=1;
	md.solidearth.settings.maxiter=10; 
	
	%time stepping:
	md.timestepping.start_time=0;
	md.timestepping.time_step=1;
	md.timestepping.final_time=1;

	%masstransport:
	md.basalforcings.groundedice_melting_rate = zeros(md.mesh.numberofvertices,1);
	md.basalforcings.floatingice_melting_rate = zeros(md.mesh.numberofvertices,1);
	md.initialization.vx = zeros(md.mesh.numberofvertices,1);
	md.initialization.vy = zeros(md.mesh.numberofvertices,1);
	md.initialization.sealevel = zeros(md.mesh.numberofvertices,1);
	md.initialization.str=0;

	md.miscellaneous.name='Transient';

	savemodel(org,md);
end % }}} 

% {{{ Static:
if perform(org,'Static'),  
	disp('   Step 4: Solve Slr solver');
	md = loadmodel(org,'Parameterization');

	md.solidearth.settings.viscous=0;
	md.solidearth.settings.selfattraction=1;
	md.solidearth.settings.elastic=1;
	md.solidearth.settings.rotation=1;

	%Physics:
	md.transient.issmb=0;
	md.transient.isstressbalance=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isslc=1;
	
	md.solidearth.requested_outputs={'Sealevel','Bed'}; 

	md=solve(md,'Transient');

	savemodel(org,md);
end % }}}

% {{{ Transient:
if perform(org,'Transient'),
	disp('   Step 5: Transient run');
	md = loadmodel(org,'Parameterization');
   	mds = loadmodel(sprintf('../icesheetmodel/ModelsSeaTest/glSlr_%d_fromobs', sigma)); %initial 1 year simulation

	[x_global, y_global]= ll2xy(md.mesh.lat, md.mesh.long, -1);

	n_steps = length(mds.results.TransientSolution);
	h0 = mds.geometry.thickness; % Initial thickness

	ice_load = zeros(md.mesh.numberofvertices, n_steps);

	disp('  Interpolating local thickness changes to global mesh...');
	for t = 1:n_steps
    	% calculate local thickness change (anomaly)
    	ht = mds.results.TransientSolution(t).Thickness;
    	dh_amundsen = ht - h0;
    	ice_load(:, t) = InterpFromMeshToMesh2d(mds.mesh.elements, mds.mesh.x,...
    	mds.mesh.y, dh_amundsen, x_global, y_global, 'default', 0);
    end
	ice_load(md.mesh.lat > 0, :) = 0;
    	
	num_time = size(ice_load,2); 
	md.masstransport.spcthickness = repmat(md.geometry.thickness,1,num_time+1); 
	md.masstransport.spcthickness(:,2:end) = md.masstransport.spcthickness(:,2:end) + ice_load;
	md.masstransport.spcthickness(end+1,:) = 0:num_time; 

	%Physics 
	md.transient.issmb=0;
	md.transient.isstressbalance=0;
	md.transient.isthermal=0;
	md.transient.ismasstransport=1;
	md.transient.isslc=1;

	md.solidearth.settings.viscous=0;
	md.solidearth.settings.selfattraction=1;
	md.solidearth.settings.elastic=1;
	md.solidearth.settings.rotation=1;
	
	%time stepping:
	md.timestepping.start_time=0;
	md.timestepping.time_step=1;
	md.timestepping.final_time=num_time; 

	md.miscellaneous.name=sprintf('glSlr_%d_THICKNESS1_Transient', sigma);
    md.settings.waitonlock = 0;
    md.verbose=verbose('solution',true,'convergence',true);
    md.cluster = mycluster;

	md.solidearth.requested_outputs = {'Sealevel','Bed'};
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);
	
	if loadonly
        %md=loadresultsfromcluster(md);
        savemodel(org,md);
    end
	
end % }}}

