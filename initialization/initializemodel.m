steps=[1:2]; loadonly=0;


codepath = '';%path to directory with ISSM binaries
execpath = '';%path to directory for simulation execution
coresrequested       = 24;
memrequested         = 5; %GB
timerequested        = 10*60; %minutes
%mycluster = pace('np',coresrequested,'login','loginname','codepath',codepath,'executionpath',execpath,'mem',memrequested,'time',timerequested);
mycluster = generic('name', oshostname(),'np',4);

org=organizer('repository','Models','prefix', 'Amundsen_','steps',steps);

%Configure ISSM
% {{{ Mesh:z
if perform(org,'Mesh'),

	disp('   -- Generating first mesh');
    md=bamg(model,'domain',['Exp/AmundsenBasinLarge.exp'],'hmax',2000);

	%Get geometry and masks
	md.geometry.bed             = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'bed');
	flag                        = find(isnan(md.geometry.bed));
	pos1                        = find(flag);
	pos2                        = find(~flag);
	md.geometry.bed(pos1)       = griddata(md.mesh.x(pos2),md.mesh.y(pos2),md.geometry.bed(pos2),md.mesh.x(pos1),md.mesh.y(pos1),'nearest');%fill in missing bed data by interpolating values from nearby points.
	md.geometry.surface         = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'surface');
	md.geometry.surface=max(md.geometry.surface,10.);
	md.geometry.base=md.geometry.bed;
	rho_water=1028.9;
	di=md.materials.rho_ice/rho_water;
	pos=find(-di/(1-di)*md.geometry.surface > md.geometry.bed);
	md.geometry.base(pos)=di/(di-1)*md.geometry.surface(pos);
	md.geometry.thickness=md.geometry.surface-md.geometry.base;
	md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;

	for i=1:1
        %Refinement loop twice
		disp('   -- Interpolating Mouginot''s velocities');
		[velx vely]=interpMouginotAnt2017(md.mesh.x,md.mesh.y);
        vel=sqrt(velx.^2+vely.^2);
        vel(find(isnan(vel)))=0;
        disp('   -- Getting thickness from MC dataset');
        thickness = interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'thickness');
		hminv=NaN(size(vel)); 
		hminv(find(thickness<20))=30000;%coarse resolution where thickness is less than 20m 
		hmaxv=NaN(size(vel));
		
		infront=ContourToMesh(md.mesh.elements,md.mesh.x,md.mesh.y,'Exp/Refinecalvingfront.exp','node',1);
		hmaxv(find(infront)) = 800;
	
		%flags=find(md.mask.ocean_levelset<=350 & md.mask.ocean_levelset>-10);
		%hmaxv(flags)=800%

       %refine mesh using surface velocities as metric
       disp('   -- remeshing with metric');
       md=bamg(md,'hmin',500,'hmax',20000,'field',vel,'err',7,'hminVertices',hminv,'hmaxVertices',hmaxv,'anisomax',3);
   end

   %Convert mesh lat long
   [md.mesh.lat,md.mesh.long]  = xy2ll(md.mesh.x,md.mesh.y,-1);
   md.mesh.epsg=3031;

   %Save model
   savemodel(org,md);
   %areas = GetAreas(md.mesh.elements,md.mesh.x,md.mesh.y);
   %plotmodel(md,'data',sqrt(2.*areas)./1000,'caxis',[1 40],'log',10);
end
% }}}

% {{{ Parameterize: 
if perform(org,'Parameterize'),
	
   %Parameterize
   md=loadmodel(org,'Mesh');
   md=setflowequation(md,'SSA','all');
   md=parameterize(md,'Par/AntarcticaAminat.par');

   %No negative thickness
   pos = find(md.geometry.thickness<=0);
   md.geometry.thickness(pos)=1; md.geometry.base(pos) = md.geometry.surface(pos)-md.geometry.thickness(pos);

   %Bed correction
   M=interpBedmachineAntarctica(md.mesh.x,md.mesh.y,'mask');
   di=md.materials.rho_ice/md.materials.rho_water;
   md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;
   pos=find(md.geometry.base<md.geometry.bed | md.mask.ocean_levelset>0);
   md.geometry.bed(pos) = md.geometry.base(pos);

	%Save OG bed
   bed = md.geometry.bed;
   bed_compare = md.geometry.bed;
   mkdir('DataFiles')
   save('DataFiles/bed.mat','bed');

	%Bed corrections
    pos = find(md.mask.ocean_levelset<0);
	md.geometry.bed(pos) = md.geometry.bed(pos)-200;
   pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/BedLowerAse.exp',2) & md.mask.ocean_levelset>-10);
   md.geometry.bed(pos) = md.geometry.bed(pos)-200;
   md.geometry.base(pos) = md.geometry.bed(pos);
   md.geometry.thickness(pos) = md.geometry.surface(pos)-md.geometry.base(pos);
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/DotsonShelf.exp',2) & md.mask.ice_levelset<0);
   md.geometry.bed(pos) = md.geometry.bed(pos)-400;
   md.geometry.base(pos) = md.geometry.bed(pos);
   md.geometry.thickness(pos) = md.geometry.surface(pos)-md.geometry.base(pos);
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/PigGlNew.exp',2) & md.mask.ocean_levelset>0);
	md.geometry.bed(pos) = md.geometry.bed(pos)-500;
   md.mask.ocean_levelset(pos) = -100;
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/ThwGlNew.exp',2) & md.mask.ocean_levelset>0);
   md.geometry.bed(pos) = md.geometry.bed(pos)-500;
   md.mask.ocean_levelset(pos) = -100;
	pos=find(ContourToNodes(md.mesh.x,md.mesh.y,'Exp/ThwGl_Upstream.exp',2) & md.mask.ocean_levelset<=0);
	md.geometry.base(pos) = md.geometry.bed(pos);
   md.geometry.thickness(pos) = md.geometry.surface(pos)-md.geometry.base(pos);
   md.mask.ocean_levelset(pos) = 100;

   %Blend edited and original beds
   bed_ed_compare = md.geometry.bed;
   bed_ed = md.geometry.bed;
   save('DataFiles/bed_ed.mat','bed_ed');
   blendbeds_tot;

	%Update ocean mask
	rho_water=1023;
	di=md.materials.rho_ice/rho_water;
	pos=find(-di/(1-di)*md.geometry.surface > md.geometry.bed);
	md.geometry.base(pos)=di/(di-1)*md.geometry.surface(pos);
	md.geometry.thickness=md.geometry.surface-md.geometry.base;
	md.mask.ocean_levelset=md.geometry.thickness+md.geometry.bed/di;

   %Save
   savemodel(org,md);

end
% }}}

% {{{ Steady-state: 
if perform(org,'Steady-state'),

	md=loadmodel(org,'Parameterize');
	mde=extrude(md,14,1.1);
	mde.timestepping.final_time=0;
	mde.timestepping.time_step=0;
	mde.thermal.isenthalpy=1;
	mde.thermal.isdynamicbasalspc=1;
	mde.thermal.stabilization=1;
	mde.initialization.watercolumn=zeros(mde.mesh.numberofvertices,1);
	mde.initialization.waterfraction=zeros(mde.mesh.numberofvertices,1);

	%First solve without advection for now
	mde.initialization.vx(:)=0;
	mde.initialization.vy(:)=0;
	mde.initialization.vz(:)=0;
	mde.initialization.vel(:)=0;

	%Set cluster and solve
	mde.miscellaneous.name='ASE_Thermal';
	mde.verbose=verbose('solution',true,'module',true,'convergence',true);
	mde.cluster = mycluster;
	mde.settings.waitonlock=0;
	mde=solve(mde,'thermal','runtimename',false,'loadonly',loadonly);
	
	if loadonly
		mde=loadresultsfromcluster(mde);
		pos=find(md.materials.rheology_B>6*10^8);
		md.materials.rheology_B(pos)=1.2*10^8;
		md.basalforcings.groundedice_melting_rate = project2d(mde,mde.results.ThermalSolution.BasalforcingsGroundediceMeltingRate,1);
		savemodel(org,md);
	end
end
% }}}
% {{{ ControlB: 
if perform(org,'ControlB'),

	%Load model
	md=loadmodel(org,'Steady-state');
   md=setflowequation(md,'SSA','all');
   md.stressbalance.restol=0.01;
   md.stressbalance.reltol=0.1;
   md.stressbalance.abstol=NaN;
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.iscontrol=1;
   md.settings.solver_residue_threshold = 1e6;

	%Extract floating elements
   mds=extract(md,md.mask.ocean_levelset<0 & md.mask.ice_levelset<0);

   %a few isolated nodes (17) are causing problems because of the extraction, so spc them
   pos=find( mds.mesh.x>-1.79*10^6 & mds.mesh.x<-1.78*10^6 & mds.mesh.y<-4.0*10^5 & mds.mesh.y>-4.2*10^5);
   mds.stressbalance.spcvx(pos)=mds.inversion.vx_obs(pos);
   mds.stressbalance.spcvy(pos)=mds.inversion.vy_obs(pos);

	%Set inversion parameters
   mds.inversion=m1qn3inversion(mds.inversion);
   mds.inversion.iscontrol=1;
   mds.inversion.maxsteps=100;
   mds.inversion.maxiter=100*10;
   mds.inversion.cost_functions=[103 101];
   mds.inversion.cost_functions_coefficients=ones(mds.mesh.numberofvertices,2);
   mds.inversion.cost_functions_coefficients(:,1)=1.5;
   mds.inversion.cost_functions_coefficients(:,2)=40;
   mds.inversion.control_parameters={'MaterialsRheologyBbar'};
   mds.inversion.min_parameters=cuffey(273.15-0.1)*ones(mds.mesh.numberofvertices,1);
   mds.inversion.max_parameters=cuffey(273.15-75)*ones(mds.mesh.numberofvertices,1);

   %Zero weight on vel_obs=0;
   pos=find(mds.inversion.vel_obs==0);
   mds.inversion.cost_functions_coefficients(pos,:)=0;

   %Solve
   mds.verbose=verbose('solution',false,'control',true);
   mds.cluster=generic('name',oshostname(),'np',24);
   mds=solve(mds,'sb');

   %Save B
   md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
   pos=find(md.materials.rheology_B>6*10^8);
   md.materials.rheology_B(pos)=1.2*10^8;
   savemodel(org,md);

end
% }}}
% {{{ Control_drag_Budd:
if perform(org,'Control_drag_Budd'),

   %Load model
   md=loadmodel(org,'ControlB');
	md.mask.ice_levelset=killberg(md);

   %Set effective pressure
   N = md.constants.g*md.materials.rho_ice*md.geometry.thickness;
   pos = find(md.geometry.base<0);
   N(pos) = max(md.constants.g*(md.materials.rho_ice*md.geometry.thickness(pos) + md.materials.rho_freshwater*md.geometry.base(pos)), 0);
   pos = find(md.mask.ocean_levelset<0 & N>0);
   N(pos) = 0;
   md.friction.effective_pressure = N;

    %Extract ice elements
   %mds=extract(md,md.mask.ice_levelset<0);

   %Set inversion controls
   md.inversion.iscontrol=1;
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.maxsteps=120;
   md.inversion.maxiter=md.inversion.maxsteps*10;
   md.inversion.cost_functions=[103 101 501];
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
   md.inversion.cost_functions_coefficients(:,1)=1;
   md.inversion.cost_functions_coefficients(:,2)=20;
   md.inversion.cost_functions_coefficients(:,3)=1*10^-8;
   pos=find(md.mask.ice_levelset>0);
   md.inversion.cost_functions_coefficients(pos,:)=0;
   %md.inversion.cost_functions_coefficients(:,3)=1*10^-6;
   pos=find(md.inversion.vel_obs<50);
   %md.inversion.cost_functions_coefficients(pos,1)=1e-2;
   md.inversion.cost_functions_coefficients(pos,2)=5;
   md.inversion.control_parameters={'FrictionCoefficient'};
   md.inversion.min_parameters=0.5*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=800*ones(md.mesh.numberofvertices,1);
   pos=find(md.mask.ocean_levelset<0);
   md.inversion.min_parameters(pos)=0;
	md.inversion.max_parameters(pos)=0;
	%Zero weight on vel_obs=1;
   pos=find(md.inversion.vel_obs==0);
   md.inversion.cost_functions_coefficients(pos,1:2)=0;
   md.transient.isgroundingline=1;
   md.groundingline.migration='SubelementMigration';
   md.verbose=verbose('control',true);
   md.cluster=generic('name',oshostname(),'np',24);

   md=solve(md,'sb');

   md.initialization.vel=(md.results.StressbalanceSolution.Vel);
   %md.friction.coefficient(md.mesh.extractedvertices)=(mds.results.StressbalanceSolution.FrictionCoefficient);
   md.friction.coefficient=md.results.StressbalanceSolution.FrictionCoefficient;
   md.friction.coefficient(find(md.mask.ocean_levelset<0))=0;
   md.inversion.iscontrol=0;

   %Save
   savemodel(org,md);

end
% }}}

% {{{ Control_drag:
if perform(org,'Control_drag')

   %Load models and set controls
   md = loadmodel(org, 'Control_drag_Budd');
   md.inversion.iscontrol = 0;
   md.transient = deactivateall(md.transient);

% {{{ Change friction coefficient to Schoof:
   %Compute friction
   N = md.constants.g*md.materials.rho_ice*md.geometry.thickness;
   pos = find(md.geometry.base<0);
   N(pos) = max(md.constants.g*(md.materials.rho_ice*md.geometry.thickness(pos) + md.materials.rho_freshwater*md.geometry.base(pos)), 0);

   %Set constants
   yts  = md.constants.yts;
   Cb   = md.friction.coefficient;
   ub   = sqrt(md.inversion.vx_obs.^2 + md.inversion.vy_obs.^2)/yts;
   taub = Cb.^2.*N.*ub;

   %Force taub and N to equal zero on floating ice
   pos = find(md.mask.ocean_levelset<0 & taub>0);
   taub(pos) = 0;
   pos = find(md.mask.ocean_levelset<0 & N>0);
   N(pos) = 0;

   %Limit taub
   taub(find(taub>3e5)) = 3e5;

	%Convert to schoof
   if 1

      %Set constants
      m = 1/3;      %prefered values: 1 or 1/n = 1/3
      Cmax = 0.8;   %prefered values between 0.17 and 0.84 ...

      %Deal with imaginary numbers
      x = min(0.999, taub./(Cmax*N));
      Cs = max(0,taub./(ub.^m .* (1 - x.^(1/m)).^m));
      Cs = sqrt(Cs);

      %Extrapolate where undefined
      flags=(Cs==0 | isnan(Cs)); pos1=find(flags); pos2=find(~flags);
      F = scatteredInterpolant(md.mesh.x(pos2),md.mesh.y(pos2),md.friction.coefficient(pos2));
      Cs(pos1) = max(0,F(md.mesh.x(pos1),md.mesh.y(pos1)));

      %Switch to Schoof
      md.friction = frictionschoof(md.friction);
      md.friction.C    = Cs;
      md.friction.Cmax = Cmax*ones(md.mesh.numberofvertices,1);
      md.friction.m    = m*ones(md.mesh.numberofelements,1);
      Cs = md.friction.C.^2;
      taub2 = Cs.*ub.^m ./ (1+ (Cs./(Cmax*max(N,1))).^(1/m) .* ub).^m;
   end

   %Solve
   md.friction.coupling = 3;
   md.friction.effective_pressure_limit = 0.05;
   md.friction.effective_pressure = N;
	%}}}
	%Extract grounded elements
   %mds=extract(md, md.mask.ice_levelset<0);
	
	%Set inversion controls
   md.inversion.iscontrol=1;
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.maxsteps=100;
   md.inversion.maxiter=md.inversion.maxsteps*10;
   md.inversion.cost_functions=[103 101 501];
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
   md.inversion.cost_functions_coefficients(:,1)=1;
   md.inversion.cost_functions_coefficients(:,2)=20;
   md.inversion.cost_functions_coefficients(:,3)=1*10^-10;
   pos=find(md.inversion.vel_obs<50);
   %md.inversion.cost_functions_coefficients(pos,1)=1e-2;
   md.inversion.cost_functions_coefficients(pos,2)=5;
   md.inversion.control_parameters={'FrictionC'};
   md.inversion.min_parameters=0.5*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=10000*ones(md.mesh.numberofvertices,1);
   pos=find(md.mask.ocean_levelset<0);
   md.inversion.min_parameters(pos)=0;
   md.inversion.max_parameters(pos)=0;
   %Zero weight on vel_obs=0;
   pos=find(md.inversion.vel_obs==0);
   md.inversion.cost_functions_coefficients(pos,1:2)=0;
   md.transient.isgroundingline=1;
   md.groundingline.migration='SubelementMigration';
   md.verbose=verbose('control',true);
   md.cluster=generic('name',oshostname(),'np',24);

   md=solve(md,'sb');

   md.friction.C=md.results.StressbalanceSolution.FrictionC;
   md.friction.C(find(md.mask.ocean_levelset<0))=0;
   md.initialization.vx=(md.results.StressbalanceSolution.Vx);
   md.initialization.vy=(md.results.StressbalanceSolution.Vy);
   md.inversion.iscontrol=0;

	%Save
   savemodel(org,md);

end
% }}}

% {{{ ControlB_helene: 
if perform(org,'ControlB_helene'),

	%Load model
	md=loadmodel(org,'Control_drag');
	
	%Redo an ice shelf inversion without extracting ice shelves
    md.inversion=m1qn3inversion(md.inversion);
    md.inversion.iscontrol=1;
	md.settings.solver_residue_threshold = 1e6;

	%Set inversion parameters
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.iscontrol=1;
   md.inversion.maxsteps=100;
   md.inversion.maxiter=100*10;
   md.inversion.cost_functions=[103 101];
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
   md.inversion.cost_functions_coefficients(:,1)=0.3;
   md.inversion.cost_functions_coefficients(:,2)=60;
	pos=find(md.mask.ice_levelset>0);
   md.inversion.cost_functions_coefficients(pos,:)=0;
   md.inversion.control_parameters={'MaterialsRheologyBbar'};
   md.inversion.min_parameters=cuffey(273.15-0.1)*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=cuffey(273.15-75)*ones(md.mesh.numberofvertices,1);
	%Do not change coefficients of grounded ice
	pos=find(md.mask.ocean_levelset>0);
	md.inversion.min_parameters(pos)=md.materials.rheology_B(pos);
	md.inversion.max_parameters(pos)=md.materials.rheology_B(pos);

   %Zero weight on vel_obs=0;
   pos=find(md.inversion.vel_obs==0);
   md.inversion.cost_functions_coefficients(pos,:)=0;

	%Solve
   md.verbose=verbose('solution',false,'control',true);
   md.cluster=generic('name',oshostname(),'np',24);
   md=solve(md,'sb');

	%Save B
   md.materials.rheology_B=md.results.StressbalanceSolution.MaterialsRheologyBbar;

	%Save results
   savemodel(org,md);

end
% }}}
% {{{ Extrusion:
if perform(org,'Extrusion')

	md=loadmodel(org,'ControlB_helene');

   %update initial velocities
   md.initialization.vx=md.results.StressbalanceSolution.Vx;
   md.initialization.vy=md.results.StressbalanceSolution.Vy;
   md.basalforcings.groundedice_melting_rate = zeros(size(md.mask.ocean_levelset));

   %Here we want to spc the velocities that are on the inflow boundary because
   % that messes up the thermal model. To use outflow however, we need the segments
   % the easiest way to do that is to extract the model, get the inflow boundary
   % and spc the velocity there.
   %pos_e = find(min(md.mask.ice_levelset(md.mesh.elements),[],2)<0);
   %flags=zeros(md.mesh.numberofvertices,1); flags(md.mesh.elements(pos_e,:))=1;
   %md2=extract(md,flags);
   %pos=find(md2.mesh.vertexonboundary & ~outflow(md2));
   %md.stressbalance.spcvx(md2.mesh.extractedvertices(pos)) = 0;
   %md.initialization.vx(  md2.mesh.extractedvertices(pos)) = 0;
   %md.stressbalance.spcvy(md2.mesh.extractedvertices(pos)) = 0;
   %md.initialization.vy(  md2.mesh.extractedvertices(pos)) = 0;
   %clear md2;

   %FIXME: don't need if parameterize now
   %md.basalforcings.floatingice_melting_rate = md.basalforcings.floatingice_melting_rate(:);

   %Extrude and set flow model as HO now
   md=extrude(md,10,1.1);
   md=setflowequation(md,'HO','all');

   savemodel(org,md);

end
% }}}
% {{{ HO:
if perform(org,'HO')

	%Load model and set inversion control
	md=loadmodel(org,'Extrusion');
    md=setflowequation(md,'HO','all');
	md.inversion.iscontrol=0;

	md.stressbalance.spcvz(:)=NaN;

	%Solve
	md.miscellaneous.name='ASE_Ho';
	md.verbose=verbose('solution',true,'module',true,'convergence',true);
	md.cluster=mycluster;
	md.settings.waitonlock=0;
	md=solve(md,'sb','runtimename',false,'loadonly',loadonly);

	%Save
	if loadonly
            %md=loadresultsfromcluster(md);
            md.initialization.vx=md.results.StressbalanceSolution.Vx;
        	md.initialization.vy=md.results.StressbalanceSolution.Vy;
        	md.initialization.vz=md.results.StressbalanceSolution.Vz;
        	md.initialization.vel=md.results.StressbalanceSolution.Vel;
        	md.initialization.pressure=md.results.StressbalanceSolution.Pressure;
                
		savemodel(org,md);
    end

end
% }}}
% {{{ Thermal:
if perform(org,'Thermal')

	%Load model and set inversion controls
    md=loadmodel(org,'HO');
    md.timestepping.time_step = 0;
    md.inversion.iscontrol = 0;
    md.thermal.isenthalpy=1;
	md.thermal.stabilization=1;

	%Change geothermal heat flux to new Stal et al. 2020 data
	md.basalforcings.geothermalflux = interpStal2020(md.mesh.x,md.mesh.y);

	%First solve without advection to get boundary condition at the divides
	md.initialization.vx(:)=0;
	md.initialization.vy(:)=0;
	md.initialization.vz(:)=0;
	md.initialization.vel(:)=0;

	%Set cluster and solve
	md.verbose=verbose('solution',false,'control',true);
    md.cluster=generic('name',oshostname(),'np',24);
	%Solve
	
	md=solve(md,'thermal');
	md=loadresultsfromcluster(md);
	md.results.Thermal0=md.results.ThermalSolution;

   %Find points on boundary and ice domain
	pos=find(md.mesh.vertexonboundary & md.mask.ice_levelset<0);
	md.thermal.spctemperature(pos)=md.results.ThermalSolution.Temperature(pos);

	md.initialization.vx=md.results.StressbalanceSolution.Vx;
	md.initialization.vy=md.results.StressbalanceSolution.Vy;
	md.initialization.vz=md.results.StressbalanceSolution.Vz;
	md.initialization.vel=md.results.StressbalanceSolution.Vel;

	md=solve(md,'thermal');
	md=loadresultsfromcluster(md);

   %Save model
   savemodel(org,md);

end
% }}}

% {{{ Control_Bnext: 
if perform(org,'Control_Bnext'),

	mde=loadmodel(org,'Thermal');
	md=collapse(mde);
	md.materials.rheology_B=cuffey(DepthAverage(mde,mde.results.ThermalSolution.Temperature));
	pos=find(md.mask.ice_levelset>0);
	md.materials.rheology_B(pos)=1.2*10^8;
	md.friction.C=mde.friction.C(1:md.mesh.numberofvertices);
	md.friction.Cmax=mde.friction.Cmax(1:md.mesh.numberofvertices);
	md.friction.m=mde.friction.m(1:md.mesh.numberofelements);
	md.friction.effective_pressure=mde.friction.effective_pressure(1:md.mesh.numberofvertices);
	md=setflowequation(md,'SSA','all');

	mds=extract(md,md.mask.ocean_levelset<0 & md.mask.ice_levelset<0.1);
	%a few isolated nodes (17) are causing problems because of the extraction again, so spc them
	pos=find( mds.mesh.x>-1.79*10^6 & mds.mesh.x<-1.78*10^6 & mds.mesh.y<-4.0*10^5 & mds.mesh.y>-4.2*10^5);
	mds.stressbalance.spcvx(pos)=mds.inversion.vx_obs(pos);
	mds.stressbalance.spcvy(pos)=mds.inversion.vy_obs(pos);
	mds.inversion=m1qn3inversion(mds.inversion);
	mds.inversion.iscontrol=1;
	mds.inversion.maxsteps=50;
	mds.inversion.maxiter=mds.inversion.maxsteps*10;
	mds.inversion.cost_functions=[103 101];
	mds.inversion.cost_functions_coefficients=ones(mds.mesh.numberofvertices,2);
	mds.inversion.cost_functions_coefficients(:,1)=1.5;
	mds.inversion.cost_functions_coefficients(:,2)=1;
	mds.inversion.control_parameters={'MaterialsRheologyBbar'};
	mds.inversion.min_parameters=paterson(273.15-0.1)*ones(mds.mesh.numberofvertices,1);
	mds.inversion.max_parameters=paterson(273.15-75)*ones(mds.mesh.numberofvertices,1);
	%Zero weight on vel_obs=0;
	pos=find(mds.inversion.vel_obs==0);
	mds.inversion.cost_functions_coefficients(pos,:)=0;

	mds.verbose=verbose('solution',false,'convergence',false);
	mds.cluster=generic('name',oshostname(),'np',24);
	mds.settings.waitonlock=1;
	mds=solve(mds,'sb');
	
	md.materials.rheology_B(mds.mesh.extractedvertices)=mds.results.StressbalanceSolution.MaterialsRheologyBbar;
	pos=find(md.mask.ice_levelset>0);
	md.materials.rheology_B(pos)=1.2*10^8;

	savemodel(org,md);
end
% }}}
% {{{ Control_dragfinal:
if perform(org,'Control_dragfinal')

   %Load models and set controls
	md=loadmodel(org,'Control_Bnext');
	md.inversion.iscontrol = 0;
    md.transient = deactivateall(md.transient);
	%Extract grounded elements
   %mds=extract(md, md.mask.ice_levelset<0);
	
	%Set inversion controls
   md.inversion.iscontrol=1;
   md.inversion=m1qn3inversion(md.inversion);
   md.inversion.maxsteps=100;
   md.inversion.maxiter=100;
   md.inversion.cost_functions=[103 101 501];
   md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,3);
   md.inversion.cost_functions_coefficients(:,1)=1;
   md.inversion.cost_functions_coefficients(:,2)=20;
   md.inversion.cost_functions_coefficients(:,3)=2*10^-10;
	pos=find(md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,:)=0;
   pos=find(md.inversion.vel_obs<50);
   md.inversion.cost_functions_coefficients(pos,2)=5;
   md.inversion.control_parameters={'FrictionC'};
   md.inversion.min_parameters=0.5*ones(md.mesh.numberofvertices,1);
   md.inversion.max_parameters=10000*ones(md.mesh.numberofvertices,1);
   pos=find(md.mask.ocean_levelset<0);
   md.inversion.min_parameters(pos)=0;
   md.inversion.max_parameters(pos)=0;
   %Zero weight on vel_obs=0;
   pos=find(md.inversion.vel_obs==0);
   md.inversion.cost_functions_coefficients(pos,1:2)=0;
	md.inversion.dfmin_frac=0.2;
   md.transient.isgroundingline=1;
   md.groundingline.migration='SubelementMigration';
   md.verbose=verbose('control',true);
   md.cluster=generic('name',oshostname(),'np',24);
   md.settings.waitonlock=1;
   md=solve(md,'sb');

   %md.friction.C(mds.mesh.extractedvertices)=(mds.results.StressbalanceSolution.FrictionC);
   md.friction.C=md.results.StressbalanceSolution.FrictionC;
   md.friction.C(find(md.mask.ocean_levelset<0))=0;
   md.inversion.iscontrol=0;

	%Save
   savemodel(org,md);

end
% }}}
% {{{ Control_Bfinal: 
if perform(org,'Control_Bfinal'),

   md = loadmodel(org, 'Control_dragfinal');

	md.inversion=m1qn3inversion(md.inversion);
	md.inversion.dfmin_frac=0.5;
	md.inversion.iscontrol=1;
	md.inversion.maxsteps=50;
	md.inversion.maxiter=50*md.inversion.maxsteps;
	md.inversion.cost_functions=[103 101];
	md.inversion.cost_functions_coefficients=ones(md.mesh.numberofvertices,2);
	md.inversion.cost_functions_coefficients(:,1)=1;
	md.inversion.cost_functions_coefficients(:,2)=60;
	md.inversion.control_parameters={'MaterialsRheologyBbar'};
	md.inversion.min_parameters=paterson(273.15-0.1)*ones(md.mesh.numberofvertices,1);
	md.inversion.max_parameters=paterson(273.15-75)*ones(md.mesh.numberofvertices,1);
	%Zero weight on vel_obs=0;
	pos=find(md.inversion.vel_obs==0 | md.mask.ice_levelset>0);
	md.inversion.cost_functions_coefficients(pos,:)=0;
	pos=find(md.mask.ocean_levelset>0);
	md.inversion.min_parameters(pos)=md.materials.rheology_B(pos);
	md.inversion.max_parameters(pos)=md.materials.rheology_B(pos);

	md.verbose=verbose('solution',false,'convergence',false);
	md.cluster=generic('name',oshostname(),'np',24);
	md=solve(md,'sb');
	md.materials.rheology_B=md.results.StressbalanceSolution.MaterialsRheologyBbar;
	pos=find(md.mask.ice_levelset>0);
	md.materials.rheology_B(pos)=1.2*10^8;

	savemodel(org,md);
end
% }}}
% {{{ Transient: 
if perform(org,'Transient')

	md=loadmodel(org,'Control_Bfinal');

	md.inversion.iscontrol=0;
	md.timestepping.start_time=0;
	md.timestepping.final_time=100;
	md.timestepping.time_step=1/12;
	md.settings.output_frequency=4;
	md.transient.isgroundingline=1;
	md.transient.isstressbalance=1;
	md.transient.ismasstransport=1;
	md.transient.isthermal=0;
	md.transient.ismovingfront=0;
	md.transient.requested_outputs={'default','IceVolume','IceVolumeAboveFloatation','GroundedArea','TotalFloatingBmb'};
	md.groundingline.migration='SubelementMigration';
	md.groundingline.friction_interpolation='SubelementFriction2';
	md.groundingline.melt_interpolation='NoMeltOnPartiallyFloating';

	%Change spcs
	md.masstransport.spcthickness=NaN*ones(md.mesh.numberofvertices,1);
	%Parameterized melting rate
	md.basalforcings=linearbasalforcings(md);
	md.basalforcings.groundedice_melting_rate=zeros(md.mesh.numberofvertices,1);
	md.basalforcings.deepwater_melting_rate=79.3397;
	md.basalforcings.deepwater_elevation=-1000;
	md.basalforcings.upperwater_elevation=-141.59;
	md.basalforcings.geothermalflux=zeros(md.mesh.numberofvertices,1);

	md.miscellaneous.name='AmundsenTransient';
	md.settings.waitonlock = 0;
	md.verbose=verbose('solution',true,'convergence',true);
	md.cluster=mycluster;
	md=solve(md,'tr','runtimename',false,'loadonly',loadonly);

	if loadonly
        %md=loadresultsfromcluster(md);
		savemodel(org,md);
	end
end
% }}}


