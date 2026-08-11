%Run coupled ice-sheet/sea-level simulation for grounding line only
%perturbation


base_dir = pwd;

F135 = fullfile(base_dir,'icesheetmodel');
F24 = fullfile(base_dir,'sealevelmodel');

% script names
S1='glcoupled_slr.m';
S2='glsealevel_icethickness1.m';
S4='glsealevel_icethickness_cycle.m';
S5='glcoupled_slr_cycle.m';

sigma = 230; %VM stress threshold

% initial ice-sheet solve from initial sea level calculated from observation

cd(F135)
steps = 1; loadonly = 0; run(S1)
steps = 1; loadonly = 1; run(S1)

% solve sea level model using initial calculated ice thickness as input
cd(F24)
steps = [1:3 5]; loadonly = 0; run(S2)
steps = 5; loadonly = 1; run(S2)

% Ice-sheet response

cd(F135)
steps = 2; loadonly = 0; run(S1)
steps = 2; loadonly = 1; run(S1)

% coupling loop

start_step = 2;
end_step = 2;

for k = start_step:end_step
    % sea-level solve
    cd(F24)
    steps = 2:4; loadonly = 0; run(S4)
    steps = 4; loadonly = 1; run(S4)


    % Ice-thickness response
    cd(F135)
    steps = 1; loadonly = 0; run(S5)
    steps = 1; loadonly = 1; run(S5)

end

fprintf('\nSimulation done\n')