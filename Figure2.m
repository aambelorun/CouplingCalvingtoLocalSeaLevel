% go to the idealized experiments output directory
cd idealized_simulations/Models

%% with zero sea level everywhere
names = {
    'cf0gl0'
    'cf5gl5'};

for i = 1:length(names)
    name = names{i};
        
    md = loadmodel(sprintf('idealized_230_%s', name));
        
    time = cell2mat({md.results.TransientSolution(:).time});
    eval(sprintf('[%s_VAF] = extract_outputs(md);', name));
    eval(sprintf('%s_sl = md.initialization.sealevel;', name));
    
    for j = 1:4
        idx = j*600;

        eval(sprintf('%s_hfront(%d).line = isoline(md,md.results.TransientSolution(%d).MaskIceLevelset,''value'',0,''output'',''matrix'');', name, j, idx));
        eval(sprintf('%s_hground(%d).line = isoline(md,md.results.TransientSolution(%d).MaskOceanLevelset,''value'',0,''output'',''matrix'');', name, j, idx));      
    end
        
    clear md
end

save('uniform_slr_outputs.mat', 'time',...
    'cf0gl0_VAF', 'cf0gl0_sl', 'cf0gl0_hfront', 'cf0gl0_hground',...
   'cf5gl5_VAF', 'cf5gl5_sl', 'cf5gl5_hfront', 'cf5gl5_hground', '-v7.3')

%%

   % 'cf1gl0'
names = {
    'cf1gl0'
    'cf5gl0'
    'cf20gl0'
    'cf0gl5'
   };

for i = 1:length(names)
    name = names{i};
    
    for k = 1:200

        if k==1
             fname = sprintf('idealized_230_%s_initial.mat', name);
        else
            fname = sprintf('idealized_230_%s_%d.mat', name, k);
        end
        
        md = loadmodel(fname);
        
        eval(sprintf('%s_t%d = cell2mat({md.results.TransientSolution(:).time});', name, k));
        eval(sprintf('[%s_VAF%d] = extract_outputs(md);', name, k));
       
        eval(sprintf('%s_hfront(%d).end = isoline(md,md.results.TransientSolution(end).MaskIceLevelset,''value'',0,''output'',''matrix'');', name, k));
        eval(sprintf('%s_hground(%d).end = isoline(md,md.results.TransientSolution(end).MaskOceanLevelset,''value'',0,''output'',''matrix'');', name, k));
        
        if k == 200
            eval(sprintf('%s_sl = md.initialization.sealevel;', name));
        end
        clear md
    end
end

save('nonuniform_slr_outputs.mat',...
    'cf1gl0_VAF*', 'cf1gl0_sl', 'cf1gl0_hfront', 'cf1gl0_hground',...
    'cf5gl0_VAF*', 'cf5gl0_sl', 'cf5gl0_hfront', 'cf5gl0_hground',...
    'cf20gl0_VAF*', 'cf20gl0_sl', 'cf20gl0_hfront', 'cf20gl0_hground',...
    'cf0gl5_VAF*', 'cf0gl5_sl', 'cf0gl5_hfront', 'cf0gl5_hground', '-v7.3')

%%
load('uniform_slr_outputs.mat')
load('nonuniform_slr_outputs.mat')
names = {
    'cf1gl0'
    'cf5gl0'
    'cf20gl0'
    'cf0gl5'
   };
%combine coupled simulations
for j = 1:length(names)
    name = names{j};

    VAFcombo = [];
    
    for kk = 1:200
        eval(sprintf('vaf = %s_VAF%d(:);', name, kk));
        VAFcombo = [VAFcombo; vaf];
    end 
    eval(sprintf('%s_VAFcombo = VAFcombo;', name));
           
end

%%
figure(2);
hold on
V0 = cf0gl0_VAF(1);
ice_density = 917;
ocean_density = 1023;

sle = ice_density/(362.5e3*ocean_density); %m SLE
time = time(1:2400);
yyaxis left
plot(time, cf0gl0_VAF(1:2400) - V0, '-', 'Color','k', LineWidth=2); 
plot(time, cf1gl0_VAFcombo - V0, '-', 'Color',[0.99 0.73 0.39], LineWidth=2); 
plot(time, cf5gl0_VAFcombo - V0, '-', 'Color',[0.95 0.45 0.10], LineWidth=2); 
plot(time, cf20gl0_VAFcombo - V0, '-', 'Color',[0.55 0.1 0.1], LineWidth=2); 
plot(time, cf0gl5_VAFcombo - V0, '-', 'Color','b', LineWidth=2); 
plot(time, cf5gl5_VAF(1:2400) - V0, '--', 'Color',[0.5 0.5 0.5],  LineWidth=2);

xlabel('Time (yrs)');
ylabel({'Ice volume above flotation (km^3)'});

grid on;
ylim([-14e4 0.5e4])
yt = [-12e4 -8e4 -4e4 0 1e4];
yticks(yt)
ax = gca;
ax.YAxis(1).Color = 'k';
set(gca, 'FontSize', 14, 'LineWidth', 1)
yl_left = ax.YAxis(1).Limits;

yyaxis right
ylim(sort(-yl_left *sle))
yticks(0:0.1:0.3)
ylabel('Sea level contribution (m SLE)')
ax.YAxis(2).Color = 'k';
ax.YAxis(2).Direction = 'reverse';

lgd = legend('No sea-level change', ...
    'Calving front only (-1 m)', ...
    'Calving front only (-5 m)', ...
    'Calving front only (-20 m)', ...
    'Grounding line only (-5 m)', ...
    'Full domain (-5 m)');

title(lgd, 'Sea-level change');
savefig(gcf, 'Figure2.fig');


%%
f=openfig('Figure2.fig');
figure(f);
exportgraphics(f, 'Figure2.eps', 'ContentType', 'vector');

%%
function [VAF] = extract_outputs(md)
    
    VAF = cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation})*10^-9;
end
