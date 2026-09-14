% compare_shorelines_BRIE.m
% Compute shoreline change = NOW.x_s - BRIE.x_s (positive = seaward if coordinates increase that way)

% clear; clc; close all;


% Files (adjust paths if needed)
fileInit = 'BRIE_initial_shoreline_Cape_Hatteras.csv';
fileNow  = 'Current_shoreline_Cape_Hatteras.csv';
outFile  = 'BRIE_shoreline_change_Cape_hatteras.csv';

% Read tables
Tinit = readtable(fileInit);
Tnow  = readtable(fileNow);

% Ensure id column is string for matching
Tinit.id = string(Tinit.id);
Tnow.id  = string(Tnow.id);

% Inner join on id to keep only common transects
T = innerjoin(Tinit, Tnow, 'Keys','id', ...
    'LeftVariables',{'id','date','x_s','UTM_E','UTM_N'}, ...
    'RightVariables',{'date','x_s','UTM_E','UTM_N'});

% Rename columns to be clear
T.Properties.VariableNames = {'id','date_init','x_init','UTM_E_init','UTM_N_init',...
    'date_now','x_now','UTM_E_now','UTM_N_now'};

% Compute change (now minus initial)
T.dx = T.x_now - T.x_init;

% Export results
writetable(T, outFile);

% Quick plots
% figure('Color','w'); 
% subplot(1,2,1);
% scatter(T.UTM_E_now, T.UTM_N_now, 40, T.dx, 'filled');
% axis equal; colorbar; colormap(jet);
% xlabel('UTM Easting'); ylabel('UTM Northing');
% title('Shoreline Change (m)');
% 
% subplot(1,2,2);
% histogram(T.dx, 30);
% xlabel('Shoreline change (m)'); ylabel('Count');
% title(sprintf('N=%d  mean=%.2f  std=%.2f', height(T), mean(T.dx,'omitnan'), std(T.dx,'omitnan')));

fprintf('Wrote %s (N=%d matching transects)\n', outFile, height(T));

%% -----------------------------
% Compare Model x_s_save(:,150) to Observed Shoreline Retreat
% -----------------------------

% Run model (ensure functions are on path)

% Ensure x_s_save numeric and has snapshot 150
x_s_save = double(b_out.x_s_save); % model stores x_s (note sign convention)
nsave = size(x_s_save,2);
if nsave < 170
    warning('Model saved only %d snapshots; using last snapshot %d instead.', nsave, nsave);
    idx_model = nsave;
else
    idx_model = 170;
end

% Convert model change to same convention as T.dx
% In barrier_model initial x_s = -BRIE_init.x_s, so BRIE x_s = -x_s_model
% Model change (now - init) in BRIE coords:
model_dx = - ( x_s_save(:,idx_model) - x_s_save(:,1) ); % nx-by-1

% Read BRIE id order used by model
BRIE_init = readtable('BRIE_initial_shoreline_Cape_Hatteras.csv');
BRIE_init.id = string(BRIE_init.id);

% Build model table and join with observed T
Tmodel = table(BRIE_init.id, model_dx, 'VariableNames', {'id','dx_model'});
T.id = string(T.id);

Tcmp = innerjoin(T, Tmodel, 'Keys','id');

% Add difference and export
Tcmp.dx_diff = Tcmp.dx_model - Tcmp.dx; % model minus observed
writetable(Tcmp, 'BRIE_shoreline_change_Cape_Hatterad.csv');

% Scatter plot model vs observed colored by point density
figure('Color','w');
x = Tcmp.dx_model;
y = Tcmp.dx;
% Compute 2D density using k-nearest neighbor (ksdensity with 2D)
xy = [x(:) y(:)];
% Use ksdensity for 2D density estimation
f = ksdensity(xy, xy);
% Sort for better plotting (low density first)
[fs, idx] = sort(f, 'ascend');
scatter(x(idx), y(idx), 36, fs, 'filled');
hold on;
lims = [min([x; y]) max([x; y])];
plot(lims, lims, 'k--','LineWidth',1);
colormap(jet); colorbar; ylabel(colorbar,'Point density');
xlabel('Model Shoreline Change (m)');
ylabel('Observed Shoreline Change (m)');
title(sprintf('Model vs Observed (snapshot %d)', idx_model));
grid on; axis equal;

fprintf('Comparison saved to BRIE_shoreline_change_with_model.csv (N=%d)\n', height(Tcmp));
%% Model vs Observed Map + Histograms
% Tcmp must contain: UTM_E_now, UTM_N_now, dx (observed), dx_model

figure('Color','w');

% Create diverging colormap red->white->blue centered at zero
ncol = 256;
top = [0 0.45 0.7];    % blue
bottom = [0.9 0.1 0.15]; % red
cmap_pos = [linspace(1,top(1),ncol/2)', linspace(1,top(2),ncol/2)', linspace(1,top(3),ncol/2)'];
cmap_neg = [linspace(bottom(1),1,ncol/2)', linspace(bottom(2),1,ncol/2)', linspace(bottom(3),1,ncol/2)'];
cmap_div = [cmap_neg; cmap_pos];
colormap(cmap_div);

% Determine symmetric color limits around zero for consistent mapping
cmax = max(abs([Tcmp.dx; Tcmp.dx_model; Tcmp.dx_diff]));
cmin = -cmax;
caxis_limits = [cmin cmax];

% Observed map
subplot(2,3,1);
scatter(Tcmp.UTM_E_now, Tcmp.UTM_N_now, 40, Tcmp.dx, 'filled');
axis equal;
caxis(caxis_limits);
hcb = colorbar; ylabel(hcb,'Shoreline change (m)');
xlabel('UTM Easting'); ylabel('UTM Northing');
title('Observed Coastsat Shoreline Change (m)');

% Model map
subplot(2,3,2);
scatter(Tcmp.UTM_E_now, Tcmp.UTM_N_now, 40, Tcmp.dx_model, 'filled');
axis equal;
caxis(caxis_limits);
hcb = colorbar; ylabel(hcb,'Shoreline change (m)');
xlabel('UTM Easting'); ylabel('UTM Northing');
title(sprintf('BRIE Shoreline Change (m)', idx_model));

% Difference map
subplot(2,3,3);
scatter(Tcmp.UTM_E_now, Tcmp.UTM_N_now, 40, Tcmp.dx_diff, 'filled');
axis equal;
caxis(caxis_limits);
hcb = colorbar; ylabel(hcb,'Shoreline change (m)');
xlabel('UTM Easting'); ylabel('UTM Northing');
title(sprintf('BRIE - Observed (m)', idx_model));

% Observed histogram with zero line
subplot(2,3,4);
h1 = histogram(Tcmp.dx, 30);
xline(0,'k--');
xlabel('Shoreline change (m)'); ylabel('Segements (n=3004)');
title(sprintf('Observed: N=%d  mean=%.2f  std=%.2f', height(Tcmp), mean(Tcmp.dx,'omitnan'), std(Tcmp.dx,'omitnan')));

% Model histogram with zero line
subplot(2,3,5);
h2 = histogram(Tcmp.dx_model, 30);
xline(0,'k--');
xlabel('Shoreline change (m)'); ylabel('Segements (n=3004)');
title(sprintf('Model: N=%d  mean=%.2f  std=%.2f', height(Tcmp), mean(Tcmp.dx_model,'omitnan'), std(Tcmp.dx_model,'omitnan')));

% Difference histogram with zero line
subplot(2,3,6);
h3 = histogram(Tcmp.dx_diff, 30);
xline(0,'k--');
xlabel('Shoreline change (m)'); ylabel('Segements (n=3004)');
title(sprintf('Difference: N=%d  mean=%.2f  std=%.2f', height(Tcmp), mean(Tcmp.dx_diff,'omitnan'), std(Tcmp.dx_diff,'omitnan')));

% Set same x-limits and y-limits for all three histograms
all_vals = [Tcmp.dx; Tcmp.dx_model; Tcmp.dx_diff];
xl = [min(all_vals,[],'omitnan') max(all_vals,[],'omitnan')];
pad = 0.05 * (xl(2)-xl(1));
xl = xl + [-pad pad];
% Determine max y across histograms by recomputing bin counts with common edges
nbins = 30;
edges = linspace(xl(1), xl(2), nbins+1);
counts1 = histcounts(Tcmp.dx, edges);
counts2 = histcounts(Tcmp.dx_model, edges);
counts3 = histcounts(Tcmp.dx_diff, edges);
ylmax = max([counts1(:); counts2(:); counts3(:)]);
for ax = [4,5,6]
    subplot(2,3,ax);
    xlim(xl);
    ylim([0 ylmax*1.05]);
end



% Optional: print basic statistics and RMSE 
% Compute residuals:
residuals = Tcmp.dx - Tcmp.dx_model; % error (observed minus model)
rmse = sqrt(mean(residuals.^2, 'omitnan'));
fprintf('Observed mean=%.2f std=%.2f | Model mean=%.2f std=%.2f | RMSE_y=%.2f (N=%d)\n', ...
    mean(Tcmp.dx,'omitnan'), std(Tcmp.dx,'omitnan'), ...
    mean(Tcmp.dx_model,'omitnan'), std(Tcmp.dx_model,'omitnan'), ...
    rmse, height(Tcmp));
