
%%

%load control run with no sea level change
mdctrl = loadmodel('idealized_simulations/Models/idealized_230_cf0gl0.mat');

tctrl = cell2mat({mdctrl.results.TransientSolution(:).time});
[ctrl_VAF] = extract_outputs(mdctrl);
ctrl_sl = mdctrl.initialization.sealevel;

for i = 1:4
    idx = i*600;
    ctrl_hfront(i).line = isoline(mdctrl,mdctrl.results.TransientSolution(idx).MaskIceLevelset,'value',0,'output','matrix');
    ctrl_hground(i).line = isoline(mdctrl,mdctrl.results.TransientSolution(idx).MaskOceanLevelset,'value',0,'output','matrix');
end

gl0=isoline(mdctrl,mdctrl.mask.ocean_levelset,'value',0,'output','matrix');%initial grounding line
cf0=isoline(mdctrl,mdctrl.mask.ice_levelset,'value',0,'output','matrix');%initial calving front

mdmesh = mdctrl.mesh;

clear mdctrl

save('control_run.mat', 'tctrl', 'ctrl_VAF', 'ctrl_sl', 'ctrl_hfront', 'ctrl_hground', 'gl0', 'cf0', 'mdmesh')

%% go to output directory
cd fullycoupled_simulations/icesheetmodel/ModelsSeaTest/

%%
 % Figure 3 

 prefixes = {'cf', 'gl', 'full'};

for i = 1:length(prefixes)
    prefix = prefixes{i};

    for k = 1:200
        if k==1
            fname = sprintf('%sSlr_230_fromobs.mat', prefix);
        else
            fname = sprintf('%sSlr_230_fromdh%d.mat', prefix, k-1);
        end 
        
        md = loadmodel(fname);
        
        eval(sprintf('%s_t%d = cell2mat({md.results.TransientSolution(:).time});', prefix, k));
        eval(sprintf('%s_h%d = cell2mat({md.results.TransientSolution(:).Thickness});', prefix, k));
        eval(sprintf('[%s_VAF%d] = extract_outputs(md);', prefix, k));
       
        eval(sprintf('%s_hfront(%d).end = isoline(md,md.results.TransientSolution(end).MaskIceLevelset,''value'',0,''output'',''matrix'');', prefix, k));
        eval(sprintf('%s_hground(%d).end = isoline(md,md.results.TransientSolution(end).MaskOceanLevelset,''value'',0,''output'',''matrix'');', prefix, k));
        
        if k == 200
            eval(sprintf('%s_sl = md.initialization.sealevel;', prefix));
        end
        clear md
    end
end

save('computed_slr_outputs.mat',...
    'cf_VAF*', 'cf_sl', 'cf_hfront', 'cf_hground',...
    'gl_VAF*', 'gl_sl', 'gl_hfront', 'gl_hground',...
    'full_VAF*', 'full_sl', 'full_hfront', 'full_hground',...
    '-v7.3')

%%
load('control_run_oldfriction.mat')
load('computed_slr_outputs.mat')
prefixes = {'cf', 'gl', 'full'};

%combine coupled simulations
for j = 1:length(prefixes)
    prefix = prefixes{j};

    VAFcombo = [];
    
    for kk = 1:200
        eval(sprintf('vaf = %s_VAF%d(:);', prefix, kk));
        VAFcombo = [VAFcombo; vaf];
    end 
    eval(sprintf('%s_VAFcombo = VAFcombo;', prefix));
           
end

%%
V0 = ctrl_VAF(1);

% plot
figure(3); clf

subplot(3, 2,[1 2])
yyaxis left
hold on
plot(tctrl, ctrl_VAF - V0, '-', 'Color', 'k', 'LineWidth', 2); 
plot(tctrl, cf_VAFcombo - V0, '-', 'Color', 'r', 'LineWidth', 2); 
plot(tctrl, gl_VAFcombo - V0, '-', 'Color', 'b', 'LineWidth', 2); 
plot(tctrl, full_VAFcombo - V0, '--', 'Color',[0.5 0.5 0.5], 'LineWidth', 2); 
grid on
xlabel('Time (yrs)');
ylabel({'Ice volume above','flotation (km^3)'});

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YTick = -6e4:2e4:0e4;
set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1)

yl_left = ax.YAxis(1).Limits;

yyaxis right

ice_density = 917;
ocean_density = 1023;
sle = ice_density/(362.5e3*ocean_density); %m SLE
ylim(sort(-yl_left * sle))
ylabel({'Sea level contribution', '(m SLE)'})
ax.YAxis(2).Color = 'k';
ax.YAxis(2).Direction = 'reverse';
legend('Control', 'Calving front','Grounding line','Full domain','Location','best');
legend boxoff
text(0.02, 0.11, '(a)', 'Units', 'normalized', 'FontSize',12, 'FontWeight','bold')

mdplot = model;
mdplot.mesh = mdmesh;
mdplot.mesh.x = mdmesh.x/1e3;
mdplot.mesh.y = mdmesh.y/1e3;
x = min(mdplot.mesh.x) + 0.02*(max(mdplot.mesh.x) - min(mdplot.mesh.x));
y = min(mdplot.mesh.y) + 0.14*(max(mdplot.mesh.y) - min(mdplot.mesh.y));

xlim_km = [-1.9e6 -1.3e6]/1e3;
ylim_km = [-750 -50];

plotmodel(mdplot, 'data', ctrl_sl, 'subplot', [3 2 3], 'caxis', [-3 0], ...
    'colorbar', 'off', 'colormap', flipud(cbrewer2('seq', 'Blues')), 'xlim', xlim_km, 'ylim', ylim_km, 'axis', 'tight equal');
ylabel('y (km)')
ax_b = gca;
text(ax_b, 0.02, 0.11, '(b)', 'Units', 'normalized', 'FontSize',12,'FontWeight','bold')

calving_cmap = [1.0 0.6 0;
        1.0 0.3 0;
        0.9 0.1 0;
        0.6 0 0];
for ii=1:4
   hf=ctrl_hfront(ii).line;
   hold on; plot(hf(:,1)/1e3,hf(:,2)/1e3,'-','Color',calving_cmap(ii,:), 'LineWidth', 2);
end

gmap = [0.6 0.6 0.6;
        0.45 0.45 0.45;
        0.3 0.3 0.3;
        0 0 0];
for ii=1:4
   hg=ctrl_hground(ii).line;
   hold on; plot(hg(:,1)/1e3,hg(:,2)/1e3,'-','Color',gmap(ii,:), 'LineWidth', 2);
end


plotmodel(mdplot, 'data', cf_sl, 'subplot', [3 2 4], 'caxis', [-3 0], ...
    'colorbar', 'off', 'colormap', flipud(cbrewer2('seq', 'Blues')),'xlim', xlim_km, 'ylim', ylim_km, 'axis', 'tight equal');
ax_c = gca;
text(ax_c, 0.02, 0.11, '(c)', 'Units', 'normalized', 'FontSize',12,'FontWeight','bold')

for ii=1:4
   kk = ii*50;
   hf=cf_hfront(kk).end;
   hold on; plot(hf(:,1)/1e3,hf(:,2)/1e3,'-','Color',calving_cmap(ii,:), 'LineWidth', 2);
end

for ii=1:4
   kk = ii*50;
   hg=cf_hground(kk).end;
   hold on; plot(hg(:,1)/1e3,hg(:,2)/1e3,'-','Color',gmap(ii,:), 'LineWidth', 2);
end

plotmodel(mdplot, 'data', gl_sl, 'subplot', [3 2 5], 'caxis', [-3 0], ...
    'colorbar', 'off', 'colormap', flipud(cbrewer2('seq', 'Blues')), 'xlim', xlim_km, 'ylim', ylim_km, 'axis', 'tight equal');
xlabel('x (km)')
ylabel('y (km)')
ax_d = gca;
text(ax_d, 0.02, 0.11, '(d)', 'Units', 'normalized', 'FontSize',12,'FontWeight','bold')

for ii=1:4
   kk = ii*50;
   hf=gl_hfront(kk).end;
   hold on; plot(hf(:,1)/1e3,hf(:,2)/1e3,'-','Color',calving_cmap(ii,:), 'LineWidth', 2);
end

for ii=1:4
   kk = ii*50;
   hg=gl_hground(kk).end;
   hold on; plot(hg(:,1)/1e3,hg(:,2)/1e3,'-','Color',gmap(ii,:), 'LineWidth', 2);
end

plotmodel(mdplot, 'data', full_sl, 'subplot', [3 2 6], 'caxis', [-3 0], ...
    'colorbar', 'off', 'colormap', flipud(cbrewer2('seq', 'Blues')), 'xlim', xlim_km, 'ylim', ylim_km, 'axis', 'tight equal');
xlabel('x (km)')
ax_e = gca;
text(ax_e, 0.02, 0.11, '(e)', 'Units', 'normalized', 'FontSize',12,'FontWeight','bold')

cb = colorbar;
cb.Label.String = 'Sea level change (m)';
clim([-3 0])
cb.Position = [0.9 0.25 0.02 0.35];
cb.FontSize = 12;
cb.Label.FontSize = 12;

%move label to left side
cb.Label.Rotation = 90;
pos = cb.Label.Position;
cb.Label.Position = [pos(1)-4.2, pos(2), pos(3)];


for ii=1:4
   kk = ii*50;
   hf=full_hfront(kk).end;
   hold on; plot(hf(:,1)/1e3,hf(:,2)/1e3,'-','Color',calving_cmap(ii,:), 'LineWidth', 2);
end

for ii=1:4
   kk = ii*50;
   hg=full_hground(kk).end;
   hold on; plot(hg(:,1)/1e3,hg(:,2)/1e3,'-','Color',gmap(ii,:), 'LineWidth', 2);
end

mapaxs = [ax_b, ax_c, ax_d, ax_e];

set(mapaxs, ...
'XTick', -1800:200:-1200, ...
'YTick', -700:200:-100, ...
'TickDir', 'in', ...
'TickLength', [0.015 0.015], ...
'Layer', 'top');

axs = findall(gcf, 'Type', 'axes');

%shift all axes to the left
for i = 1:length(axs)
    axs(i).FontSize = 12;
    pos = axs(i).Position;
    pos(1) = pos(1) - 0.02;
    axs(i).Position = pos; 
end

% add legend to plot b
ax_map = axes('Position', [0.18 0.015 0.65 0.04],'Visible','off');
hold(ax_map, 'on')

h_cf_label =  plot(ax_map, nan, nan, 'LineStyle', 'none');
h_gl_label =  plot(ax_map, nan, nan, 'LineStyle', 'none');

h_cf = gobjects(4,1);
h_gl = gobjects(4,1);

for ii=1:4
    h_cf(ii) = plot(ax_map, nan, nan, '-', 'Color', calving_cmap(ii,:), 'LineWidth', 2);
    
    h_gl(ii) = plot(ax_map, nan, nan, '-', 'Color', gmap(ii,:), 'LineWidth', 2);
end
    
leg = legend(ax_map, [h_cf_label; h_cf; h_gl_label; h_gl], ...
    'Calving front position', '50 yr','100 yr','150 yr','200 yr', ...
    'Grounding line position', '50 yr','100 yr','150 yr','200 yr', ...
    'Orientation', 'horizontal', ...
    'NumColumns', 5, ...
    'Location','northeast');
legend boxoff

leg.FontSize = 12;
savefig(gcf, 'Figure3_newmask.fig');


%%
f=openfig('Figure3.fig');
figure(f);
exportgraphics(f, 'Figure3.eps', 'ContentType', 'vector');

%%
function [VAF] = extract_outputs(md)
    
    VAF = cell2mat({md.results.TransientSolution(:).IceVolumeAboveFloatation})*10^-9;
end
