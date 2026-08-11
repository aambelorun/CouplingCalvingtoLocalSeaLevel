% Run idealized SLR experiment locally
% Equivalent to cluster shell script

F135 = pwd;

% script name
S1='idealizedslr_cycle.m';
	
% set sea level perturbation value
cf = 0; %calving front perturbation
gl = 0; %grounding line perturbation

sigma = 230; %set VM stress threshold

cd(F135)

fprintf('Starting simulation\n')

if cf == gl
   steps = 3;loadonly = 0;
   run(S1)

   %load result
   steps = 3;loadonly = 1;
   run(S1)

   % Restart annually for cf and gl perturbation simulations
else
	steps = 1; loadonly = 0;
    run(S1)

    steps = 1; loadonly = 1;
    run(S1)

    %restarting annually loop

	start_step = 2;
    end_step = 3;

	for k = start_step:end_step

        %Restart from previous time step
        steps = 2; loadonly = 0;
        run(S1)

        steps = 2; loadonly = 1;
        run(S1)
    end

end

fprintf('\nsimulation done\n')
