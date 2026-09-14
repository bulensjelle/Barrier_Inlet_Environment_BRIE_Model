function write_mipNC_submission(brie_output_file, brie_initial_shoreline_file, ...
    geojson_file, template_csv_file, output_csv_file, id_start, id_end, model_start_year)
% WRITE_MIPNC_SUBMISSION_COASTSAT  Invert BRIE shoreline output (x_s_save)
% back to each CoastSat transect's own shoreline-distance convention, and
% write the mipNC submission CSV, values on 1 January of each year only.
%
% This mirrors (and inverts) the geometry from the CoastSat->BRIE
% conversion script: BRIE_xs = c_i + m_i*xs0 per transect i, where c_i
% and m_i depend only on transect/baseline geometry. So:
%     xs0 = (BRIE_xs - c_i) / m_i
%
% USAGE
%   write_mipNC_submission_coastsat('b_out.mat', ...
%       'BRIE_initial_shoreline.csv', 'NC_transects.geojson', ...
%       'mipNC_1871-2012_SubmissionNameGoesHere.csv', ...
%       'mipNC_1871-2012_MySubmission.csv', ...
%       "usa_NC_0049_0230", "usa_NC_0032_0001", 1850)
%
% INPUTS
%   brie_output_file           - .mat with b_out struct (or bare
%                                 variables) containing x_s_save
%                                 [ny x n_saved_years], meters.
%   brie_initial_shoreline_file- the CSV your model reads at init
%                                 (has an 'id' column in the SAME order
%                                 as the model's alongshore grid, rows
%                                 1:ny -> grid index 1:ny). This is what
%                                 lets us match template columns to grid
%                                 rows by transect ID.
%   geojson_file                - the CoastSat transects geojson used in
%                                 your conversion script (gives each
%                                 transect's endpoints, in lon/lat).
%   template_csv_file           - submission template (column order +
%                                 full daily date range).
%   output_csv_file             - path to write the filled submission.
%   id_start, id_end             - the two transect IDs that define the
%                                 BRIE baseline (must match what you used
%                                 in the CoastSat->BRIE script). Pass as
%                                 string, e.g. "usa_NC_0049_0230".
%   model_start_year             - calendar year of x_s_save column 1.
%                                 Defaults to 1871.
%
% NOTES / ASSUMPTIONS
%   - UTM zone is hardcoded to EPSG:32618 (Cape Hatteras / UTM 18N), same
%     as your conversion script. Change utmCRS below if a different run
%     uses a different zone.
%   - Transect origin/direction are taken from the FIRST and LAST
%     geojson coordinate of each transect, same as the forward script.
%   - Output values are xs0 in each transect's own native distance
%     convention (same units/meaning as the historical CoastSat CSVs'
%     "Var2" column) - NOT UTM coordinates. Let me know if you also want
%     UTM_E/UTM_N per point; that needs a wider (2-col-per-transect)
%     format so I've left it out of this single-value-per-cell CSV.
%   - Transects that run nearly parallel to the baseline (m_i close to
%     0) make the inversion unstable - these are flagged as warnings and
%     left blank in the output rather than producing a garbage number.

if nargin < 8 || isempty(model_start_year)
    model_start_year = 1871;
end
id_start = string(id_start);
id_end   = string(id_end);

M_MIN = 0.05; % below this |m_i|, treat inversion as unreliable (near-parallel transect)

%% 1. Load BRIE model output
S = load(brie_output_file);
if isfield(S, 'b_out')
    b_out = S.b_out;
else
    b_out = S;
end
if ~isfield(b_out, 'x_s_save')
    error('Field "x_s_save" not found in %s.', brie_output_file);
end
x_s_save = double(b_out.x_s_save); % [ny x n_saved_years]
ny = size(x_s_save, 1);
n_years_saved = size(x_s_save, 2);
model_years = model_start_year : (model_start_year + n_years_saved - 1);

%% 2. Load BRIE_initial_shoreline.csv for the grid-index <-> transect ID order
BRIE_init = readtable(brie_initial_shoreline_file);
BRIE_init.id = string(BRIE_init.id);
if height(BRIE_init) ~= ny
    warning(['BRIE_initial_shoreline has %d rows but x_s_save has ny=%d. ' ...
             'They should match 1:1 by row (grid index) - double check ' ...
             'this is the same run/config.'], height(BRIE_init), ny);
end
grid_id = BRIE_init.id; % grid_id(i) = transect ID at alongshore grid index i

%% 3. Load geojson and build per-transect origin + unit direction (UTM)
G = jsondecode(fileread(geojson_file));
n_all = length(G.features);
ID_all = strings(n_all,1);
LAT_all = cell(n_all,1);
LON_all = cell(n_all,1);
for i = 1:n_all
    ID_all(i) = string(G.features(i).properties.id);
    coords = G.features(i).geometry.coordinates;
    LON_all{i} = coords(:,1);
    LAT_all{i} = coords(:,2);
end

utmCRS = projcrs(32618); % EPSG:32618, Cape Hatteras UTM 18N - change if needed

tx1_all = nan(n_all,1); ty1_all = nan(n_all,1);
dx_all  = nan(n_all,1); dy_all  = nan(n_all,1);
for i = 1:n_all
    [x,y] = projfwd(utmCRS, LAT_all{i}, LON_all{i});
    tx1_all(i) = x(1);
    ty1_all(i) = y(1);
    ddx = x(end)-x(1);
    ddy = y(end)-y(1);
    dL = sqrt(ddx^2+ddy^2);
    dx_all(i) = ddx/dL;
    dy_all(i) = ddy/dL;
end

%% 4. Recompute the BRIE baseline (must match the forward conversion exactly)
i1 = find(ID_all == id_start, 1);
i2 = find(ID_all == id_end, 1);
if isempty(i1) || isempty(i2)
    error('Baseline transect IDs not found in geojson.');
end
x1 = tx1_all(i1); y1 = ty1_all(i1);
x2 = tx1_all(i2); y2 = ty1_all(i2);
bx = x2-x1; by = y2-y1;
L = sqrt(bx^2+by^2);
ux = bx/L; uy = by/L;   % alongshore unit vector
nx = -uy;  nyv = ux;    % cross-shore (normal) unit vector

%% 5. Per-transect c_i, m_i (indexed to match grid_id order, i.e. grid index)
c_i = nan(ny,1);
m_i = nan(ny,1);
for gi = 1:ny
    match = find(ID_all == grid_id(gi), 1);
    if isempty(match)
        warning('Transect %s (grid index %d) not found in geojson - will be left blank.', grid_id(gi), gi);
        continue
    end
    c_i(gi) = (tx1_all(match)-x1)*nx + (ty1_all(match)-y1)*nyv;
    m_i(gi) = dx_all(match)*nx + dy_all(match)*nyv;
end

unstable = abs(m_i) < M_MIN;
if any(unstable)
    warning(['%d transect(s) are nearly parallel to the baseline (|m_i| < %.2f) - ' ...
             'inversion is unreliable for these and they will be left blank.'], ...
             sum(unstable), M_MIN);
end

%% 6. Invert every saved model year: xs0 = (BRIE_xs - c_i) / m_i
xs0_save = nan(ny, n_years_saved);
valid = ~isnan(c_i) & ~isnan(m_i) & ~unstable;
xs0_save(valid, :) = (-x_s_save(valid, :) - c_i(valid)) ./ m_i(valid);  

%% 7. Read template: column order (transect IDs) and full date range
fid = fopen(template_csv_file, 'r');
header_line = fgetl(fid);
template_cols = strsplit(header_line, ',');
transect_ids = string(template_cols(2:end));
n_transects = numel(transect_ids);

all_dates_str = {};
line = fgetl(fid);
while ischar(line)
    c = strsplit(line, ',');
    all_dates_str{end+1} = c{1}; %#ok<AGROW>
    line = fgetl(fid);
end
fclose(fid);
all_dates = datetime(all_dates_str, 'InputFormat', 'yyyy-MM-dd');
n_days = numel(all_dates);

%% 8. Match template columns to grid rows BY ID (not by position)
col_to_grid = nan(n_transects,1);
for c = 1:n_transects
    match = find(grid_id == transect_ids(c), 1);
    if isempty(match)
        warning('Template transect %s not found in BRIE_initial_shoreline - left blank.', transect_ids(c));
        continue
    end
    col_to_grid(c) = match;
end

%% 9. Assemble output: NaN everywhere except Jan-1 rows with model data
out = nan(n_days, n_transects);
is_jan1 = (month(all_dates) == 1) & (day(all_dates) == 1);
jan1_years = year(all_dates(is_jan1));
jan1_row_idx = find(is_jan1);

for k = 1:numel(jan1_row_idx)
    yr = jan1_years(k);
    yr_col = find(model_years == yr, 1);
    if isempty(yr_col)
        continue
    end
    row = jan1_row_idx(k);
    good = ~isnan(col_to_grid);
    out(row, good) = xs0_save(col_to_grid(good), yr_col)';
end

%% 10. Write CSV, reusing the template header verbatim
date_strs = cellstr(datestr(all_dates, 'yyyy-mm-dd'));
fid = fopen(output_csv_file, 'w');
fprintf(fid, '%s\n', header_line);
for r = 1:n_days
    vals = out(r, :);
    strvals = cell(1, n_transects);
    for c = 1:n_transects
        if isnan(vals(c))
            strvals{c} = '';
        else
            strvals{c} = sprintf('%.3f', vals(c));
        end
    end
    fprintf(fid, '%s,%s\n', date_strs{r}, strjoin(strvals, ','));
end
fclose(fid);

fprintf('Wrote %s (%d rows, %d transects, %d populated Jan-1 rows)\n', ...
    output_csv_file, n_days, n_transects, sum(is_jan1 & any(~isnan(out),2)));

end

write_mipNC_submission_coastsat('b_out.mat', ...
     'BRIE_initial_shoreline_Cape_Hatteras.csv', 'NC_transects.geojson', ...
     'C:\Users\6679242\Documents\Shoreshop3\Submissions\mipNC_1871-2012_BRIE.csv', ... ,
     'C:\Users\6679242\Documents\Shoreshop3\Submissions\mipNC_1871-2012_BRIE_Cape_Hatteras.csv', ...     
     "usa_NC_0049_0230", "usa_NC_0032_0001", 1850)
