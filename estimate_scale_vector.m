function [scaleVec, diagnostics] = estimate_scale_vector( ...
    optVars, paramNames, nSamples, nReplicatesPerSample)
%ESTIMATE_SCALE_VECTOR Estimate objective-function scales from LHS samples.
%
% Example:
%   [scaleVec, diagnostics] = estimate_scale_vector(optVars, paramNames, 30, 1);
%
% Inputs
%   optVars               Optimisation-variable array with .Range and .Transform
%   paramNames            Cell array of matching b_struct field names
%   nSamples              Number of Latin-hypercube samples (default = 30)
%   nReplicatesPerSample  Number of stochastic model repeats (default = 1)
%
% Outputs
%   scaleVec              Scale for:
%                         [Nmean_diff, Nmed_diff, Nrmse, ks_stat, Nwasserstein]
%   diagnostics           Structure containing samples, metrics, and failure info
%
% IMPORTANT:
% Add the following to barrier_model_buffer.m:
%
%   drowned = false;
%
% before the main time loop, and:
%
%   b_out.drowned = drowned;
%
% near the end of the function, before its final "end".

%% Defaults and checks
if nargin < 3 || isempty(nSamples)
    nSamples = 30;
end

if nargin < 4 || isempty(nReplicatesPerSample)
    nReplicatesPerSample = 1;
end

paramNames = cellstr(paramNames);
nParams = numel(paramNames);

if numel(optVars) ~= nParams
    error('optVars and paramNames must have the same number of elements.');
end

metricNames = {'Nmean_diff','Nmed_diff','Nrmse','ks_stat','Nwasserstein'};
nMetrics = numel(metricNames);

%% Get parameter bounds and transformations
lb = zeros(1,nParams);
ub = zeros(1,nParams);
isLog = false(1,nParams);

for k = 1:nParams
    v = optVars(k);

    lb(k) = v.Range(1);
    ub(k) = v.Range(2);
    isLog(k) = strcmpi(v.Transform,'log');

    if ~isfinite(lb(k)) || ~isfinite(ub(k)) || lb(k) >= ub(k)
        error('Invalid bounds for parameter %s.', paramNames{k});
    end

    if isLog(k) && lb(k) <= 0
        error('Log-transformed parameter %s must have a positive lower bound.', ...
            paramNames{k});
    end
end

%% Latin-hypercube parameter sampling
lhs = lhsdesign(nSamples,nParams);
X = zeros(nSamples,nParams);

for k = 1:nParams
    if isLog(k)
        X(:,k) = 10.^(log10(lb(k)) + lhs(:,k) .* ...
            (log10(ub(k)) - log10(lb(k))));
    else
        X(:,k) = lb(k) + lhs(:,k) .* (ub(k) - lb(k));
    end
end

%% Allocate outputs
allMetrics = nan(nSamples,nMetrics);

nValidReps   = zeros(nSamples,1);
nDrownedReps = zeros(nSamples,1);
nErrorReps   = zeros(nSamples,1);

messages = cell(nSamples,1);

%% Run Latin-hypercube samples
for i = 1:nSamples

    % Fresh parameter structure for every sample
    b_struct = initialize_barrier_model_Cape_Hatteras();

    % Apply sampled parameter values
    for k = 1:nParams
        if ~isfield(b_struct,paramNames{k})
            error('Parameter "%s" is not a field in b_struct.', paramNames{k});
        end

        b_struct.(paramNames{k}) = X(i,k);
    end

    reps = nan(nReplicatesPerSample,nMetrics);

    for r = 1:nReplicatesPerSample
        try
            b_out = barrier_model_buffer(b_struct);

            % A drowned barrier is not a successful model simulation.
            if ~isfield(b_out,'drowned')
                error(['b_out.drowned is missing. Add "drowned = false" before ', ...
                    'the time loop and "b_out.drowned = drowned" near the ', ...
                    'end of barrier_model_buffer.']);
            end

            if b_out.drowned
                nDrownedReps(i) = nDrownedReps(i) + 1;

                messages{i} = sprintf('%s Rep %d drowned.', ...
                    messages{i},r);

                continue
            end

            % Retrieve calibration metrics
            metrics = nan(1,nMetrics);

            for m = 1:nMetrics
                if ~isfield(b_out,metricNames{m})
                    error('b_out.%s does not exist.',metricNames{m});
                end

                metrics(m) = b_out.(metricNames{m});
            end

            % Reject non-finite output metrics
            if any(~isfinite(metrics))
                nErrorReps(i) = nErrorReps(i) + 1;

                messages{i} = sprintf('%s Rep %d returned NaN/Inf metrics.', ...
                    messages{i},r);

                continue
            end

            reps(r,:) = metrics;
            nValidReps(i) = nValidReps(i) + 1;

        catch ME
            nErrorReps(i) = nErrorReps(i) + 1;

            messages{i} = sprintf('%s Rep %d error: %s', ...
                messages{i},r,ME.message);

            warning('Sample %d, replicate %d failed: %s',i,r,ME.message);
        end
    end

    % Conservatively retain a sample only if ALL replicates succeeded.
    % This matters when nReplicatesPerSample > 1.
    validReps = all(isfinite(reps),2);

    if all(validReps) && nDrownedReps(i) == 0 && nErrorReps(i) == 0
        allMetrics(i,:) = mean(reps,1);

        fprintf('Sample %d/%d valid (%d replicate(s)).\n', ...
            i,nSamples,nReplicatesPerSample);
    else
        fprintf('Sample %d/%d invalid: %d valid, %d drowned, %d error.\n', ...
            i,nSamples,nValidReps(i),nDrownedReps(i),nErrorReps(i));
    end
end

%% Retain only successful parameter sets
validSamples = all(isfinite(allMetrics),2) & ...
               nDrownedReps == 0 & ...
               nErrorReps == 0;

metricsOK = allMetrics(validSamples,:);

fprintf('\n--- Scale-vector summary ---\n');
fprintf('Valid samples:   %d / %d\n',sum(validSamples),nSamples);
fprintf('Drowned samples: %d / %d\n',sum(nDrownedReps > 0),nSamples);
fprintf('Error samples:   %d / %d\n',sum(nErrorReps > 0),nSamples);

if size(metricsOK,1) < 2
    error(['Fewer than two valid simulations remain. Increase nSamples or ', ...
        'narrow the parameter ranges to avoid widespread drowning.']);
end

if size(metricsOK,1) < 10
    warning(['Only %d valid samples remain. This is sufficient for debugging, ', ...
        'but not for a stable calibration scale estimate.'],size(metricsOK,1));
end

%% Calculate scales
% MAD is robust to occasional extreme valid simulations.
madScale = 1.4826 .* mad(metricsOK,1,1);

% Used only if MAD = 0 but valid simulations have nonzero variability.
stdScale = std(metricsOK,0,1);

% Useful diagnostic, not directly used unless inspecting results.
rangeScale = max(metricsOK,[],1) - min(metricsOK,[],1);

scaleVec = madScale;

% If MAD is zero because many values are duplicated, use standard deviation.
useStdFallback = (scaleVec == 0 | ~isfinite(scaleVec)) & ...
                 isfinite(stdScale) & stdScale > 0;

scaleVec(useStdFallback) = stdScale(useStdFallback);

% If a metric is truly constant, assign scale = 1 to prevent divide-by-zero.
% Such a metric cannot discriminate among parameter sets.
constantMetric = (scaleVec == 0 | ~isfinite(scaleVec));
scaleVec(constantMetric) = 1;

%% Print diagnostics
fprintf('\nMetric scaling diagnostics:\n');

for m = 1:nMetrics
    fprintf(['%-14s | MAD = %12.6g | std = %12.6g | range = %12.6g ', ...
             '| final scale = %12.6g'], ...
             metricNames{m},madScale(m),stdScale(m), ...
             rangeScale(m),scaleVec(m));

    if useStdFallback(m)
        fprintf('  [std fallback]');
    end

    if constantMetric(m)
        fprintf('  [constant metric: scale set to 1]');
    end

    fprintf('\n');
end

fprintf('\nscaleVec = [');
fprintf(' %.8g',scaleVec);
fprintf(' ]\n');

%% Return diagnostics
sampleTable = array2table(X,'VariableNames',paramNames);

for m = 1:nMetrics
    sampleTable.(metricNames{m}) = allMetrics(:,m);
end

sampleTable.validSample = validSamples;
sampleTable.nValidReps = nValidReps;
sampleTable.nDrownedReps = nDrownedReps;
sampleTable.nErrorReps = nErrorReps;
sampleTable.message = messages;

diagnostics = struct();
diagnostics.metricNames = metricNames;
diagnostics.parameterSamples = X;
diagnostics.allMetrics = allMetrics;
diagnostics.validSamples = validSamples;
diagnostics.metricsOK = metricsOK;
diagnostics.madScale = madScale;
diagnostics.stdScale = stdScale;
diagnostics.rangeScale = rangeScale;
diagnostics.usedStdFallback = useStdFallback;
diagnostics.constantMetric = constantMetric;
diagnostics.sampleTable = sampleTable;

end