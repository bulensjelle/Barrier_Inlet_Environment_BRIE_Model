%% ============================================================
% CoastSat -> BRIE Initial Shoreline
%
% Cape Hatteras
% UTM Zone 18N (EPSG:32618)
%
% Baseline between two transects = xs 0
%
% Shoreline value taken from the CoastSat date CLOSEST to a
% target reference date (default: 1 Jan 2024)
%
% ============================================================

clear; clc; close all;


%% -----------------------------
% INPUTS
% ------------------------------

geojsonFile = ...
'C:\Users\6679242\Documents\Shoreshop3\Data\NC_transects.geojson';


csvFolder = ...
'C:\Users\6679242\Documents\Shoreshop3\Data\NC_Coastsat_smoothed';


% Choose two CoastSat transects defining BRIE baseline

id_start = "usa_NC_0049_0230";
id_end   = "usa_NC_0032_0001";


% Reference date: for each transect, the CoastSat sample whose
% date is closest to this value is used as the "initial" shoreline

targetDate = datetime(2024,1,1);


outputFile = ...
'BRIE_2024_shoreline_Cape_hatteras.csv';

%% -----------------------------
% READ GEOJSON
% -----------------------------

G = jsondecode(fileread(geojsonFile));

n_all = length(G.features);

ID_all  = strings(n_all,1);
LAT_all = cell(n_all,1);
LON_all = cell(n_all,1);


for i = 1:n_all

    % ID from GeoJSON
    ID_all(i) = string(G.features(i).properties.id);


    % Coordinates
    coords = G.features(i).geometry.coordinates;

    % GeoJSON order = longitude, latitude
    LON_all{i} = coords(:,1);
    LAT_all{i} = coords(:,2);

end



%% -----------------------------
% SELECT ONLY usa_NC_0049_xxxx TO usa_NC_0032_xxxx
% -----------------------------

segment_number = zeros(n_all,1);


for i = 1:n_all

    parts = split(ID_all(i),"_");

    % ID format:
    % usa_NC_0049_0001
    %
    % parts{1}=usa
    % parts{2}=NC
    % parts{3}=0049
    % parts{4}=0001

    segment_number(i) = str2double(parts{3});

end


% Keep only segments 32-49

idx_keep = segment_number >= 32 & segment_number <= 49;


% Apply selection

ID  = ID_all(idx_keep);
LAT = LAT_all(idx_keep);
LON = LON_all(idx_keep);


n = length(ID);



% Optional: sort from 0049 -> 0032

[~,order] = sort(segment_number(idx_keep),'descend');

ID  = ID(order);
LAT = LAT(order);
LON = LON(order);


disp("Selected transects:")
disp(ID)

%% -----------------------------
% ORDER TRANSECTS FOR BRIE
% 0049 -> 0032
% XXXX high -> low within each group
% -----------------------------

segment_group = zeros(length(ID),1);
segment_number = zeros(length(ID),1);


for i = 1:length(ID)

    parts = split(ID(i),"_");

    % Example:
    % usa_NC_0049_0123
    %
    % parts{3} = 0049
    % parts{4} = 0123

    segment_group(i) = str2double(parts{3});
    segment_number(i) = str2double(parts{4});

end


% Sort:
% group: descending
% segment: descending

sort_table = table(segment_group,...
                   segment_number,...
                   (1:length(ID))');


sort_table = sortrows(sort_table,...
    {'segment_group','segment_number'},...
    {'descend','descend'});


idx_order = sort_table.Var3;


% Apply ordering

ID  = ID(idx_order);
LAT = LAT(idx_order);
LON = LON(idx_order);


n = length(ID);


disp("BRIE order:")
disp(ID)

%% -----------------------------
% CONVERT TRANSECTS TO UTM
% -----------------------------

% Cape Hatteras
utmCRS = projcrs(32618);


UTM_X = cell(n,1);
UTM_Y = cell(n,1);



for i = 1:n


    [x,y] = projfwd(utmCRS,...
                    LAT{i},...
                    LON{i});


    UTM_X{i}=x;
    UTM_Y{i}=y;


end



%% -----------------------------
% CREATE BRIE BASELINE
% -----------------------------

i1 = find(ID==id_start);
i2 = find(ID==id_end);


if isempty(i1) || isempty(i2)

    error('Baseline IDs not found')

end



% Baseline points

x1 = UTM_X{i1}(1);
y1 = UTM_Y{i1}(1);


x2 = UTM_X{i2}(1);
y2 = UTM_Y{i2}(1);



% Baseline vector

bx = x2-x1;
by = y2-y1;


L = sqrt(bx^2+by^2);



% Alongshore unit vector

ux = bx/L;
uy = by/L;



% Cross-shore unit vector

nx = -uy;
ny = ux;



%% -----------------------------
% COASTSAT SHORELINE EXTRACTION
% (value closest to targetDate)
% -----------------------------


BRIE_xs = nan(n,1);

shore_X = nan(n,1);
shore_Y = nan(n,1);

initial_date = NaT(n,1);

date_offset_days = nan(n,1);   % how far the chosen sample is from targetDate



for i = 1:n


    id = ID(i);


    csvFile = fullfile(csvFolder,id+".csv");


    if ~isfile(csvFile)

        warning("Missing CSV: %s",id)
        continue

    end



    % Read CSV without headers

    C = readtable(csvFile,...
        'ReadVariableNames',false);



    dates = datetime(C.Var1);

    xs = C.Var2;



    % Sort oldest first (not required for nearest-date search,
    % but kept for consistency / easier debugging)

    [dates,idx] = sort(dates);

    xs = xs(idx);



    % Find the sample closest to targetDate

    [min_diff,idx_closest] = min(abs(dates - targetDate));

    xs0 = xs(idx_closest);


    initial_date(i)     = dates(idx_closest);
    date_offset_days(i) = days(min_diff);



    %% Convert CoastSat distance to UTM


    tx = UTM_X{i};
    ty = UTM_Y{i};



    % Transect direction

    dx = tx(end)-tx(1);
    dy = ty(end)-ty(1);


    tL = sqrt(dx^2+dy^2);


    dx = dx/tL;
    dy = dy/tL;



    % Shoreline UTM position

    sx = tx(1)+xs0*dx;
    sy = ty(1)+xs0*dy;



    shore_X(i)=sx;
    shore_Y(i)=sy;



    %% Convert to BRIE xs coordinate


    vx = sx-x1;
    vy = sy-y1;


    % perpendicular distance from baseline

    BRIE_xs(i)=vx*nx + vy*ny;



end


% Warn about any transects whose closest sample is far from targetDate
far_idx = date_offset_days > 90;   % flag if >90 days away, adjust as needed
if any(far_idx)
    warning("%d transect(s) have no CoastSat sample within 90 days of target date:",...
        sum(far_idx))
    disp(ID(far_idx))
end



%% -----------------------------
% EXPORT
% -----------------------------


BRIE = table(...
    ID,...
    initial_date,...
    date_offset_days,...
    BRIE_xs,...
    shore_X,...
    shore_Y,...
    ...
    'VariableNames',...
    {'id','date','days_from_target','x_s','UTM_E','UTM_N'});


writetable(BRIE,outputFile)



disp("Finished")
disp(outputFile)



%% -----------------------------
% PLOT CHECK
% -----------------------------


figure('Color','w')
hold on


% Transects

for i=1:n

    plot(UTM_X{i},UTM_Y{i},...
        'Color',[0.7 0.7 0.7])

end



% BRIE baseline

plot([x1 x2],...
     [y1 y2],...
     'r-',...
     'LineWidth',3)



% Initial shoreline

scatter(shore_X,...
        shore_Y,...
        30,...
        BRIE_xs,...
        'filled')


axis equal
colorbar


xlabel('UTM Easting')
ylabel('UTM Northing')


title('CoastSat shoreline closest to 1 Jan 2024, in BRIE coordinates')


legend('CoastSat transects',...
       'BRIE baseline',...
       'Shoreline @ target date')