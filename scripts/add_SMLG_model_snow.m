function add_SMLG_model_snow(repoRoot)

close all

finalDir = fullfile(repoRoot, 'data', 'final');

era5Dir = fullfile(repoRoot, ...
    'data', 'model', 'ERA5');
merraDir = fullfile(repoRoot, ...
    'data', 'model', 'MERRA2');

inputFile = fullfile(finalDir, ...
    'snow_model_matchups.mat');

outputFile = fullfile(finalDir, ...
    'snow_model_matchups_with_models.mat');

ncfile_era5 = fullfile(era5Dir, ...
    'SM_snod_ERA5_01Aug1980-31Jul2021_v01.nc');

ncfiles_merra = { ...
    fullfile(merraDir,'SM_snod_MERRA2_ease_01Aug1980-31Jul2018.nc')
    fullfile(merraDir,'SM_snod_MERRA2_ease_01Aug2018-31Jul2021.nc')
    fullfile(merraDir,'SM_snod_MERRA2_ease_01Aug2021-31Jul2022.nc')
    fullfile(merraDir,'SM_snod_MERRA2_ease_01Aug2022-31Jul2023.nc')
    };

cutoff_date_era5  = datetime(2021,7,31);
cutoff_date_merra = datetime(2023,7,31);
n_clim_years = 10;

time_half = 2;   % +/- 2 days
grid_half = 2;   % 5x5 cells

load(inputFile, 'D');

if ~exist('D','var') || ~istable(D)
    error('Variable D was not found in %s or is not a table.', inputFile);
end

required_vars = ["date","lat","lon","hs"];
missing_vars = required_vars(~ismember(required_vars, string(D.Properties.VariableNames)));
if ~isempty(missing_vars)
    error('Table D is missing required variables: %s', strjoin(cellstr(missing_vars), ', '));
end

D.date = ensure_datetime_column(D.date, height(D));
D.lat  = ensure_numeric_column(D.lat,  height(D));
D.lon  = ensure_numeric_column(D.lon,  height(D));
D.hs   = ensure_numeric_column(D.hs,   height(D));
D.lon  = wrapTo180_local(D.lon);

if ~ismember("hs_smlg", string(D.Properties.VariableNames))
    D.hs_smlg = nan(height(D),1);
else
    D.hs_smlg = ensure_numeric_column(D.hs_smlg, height(D));
    D.hs_smlg(:) = NaN;
end

if ~ismember("t_smlg", string(D.Properties.VariableNames))
    D.t_smlg = NaT(height(D),1);
else
    D.t_smlg = ensure_datetime_column(D.t_smlg, height(D));
    D.t_smlg(:) = NaT;
end

if ~ismember("hs_merra", string(D.Properties.VariableNames))
    D.hs_merra = nan(height(D),1);
else
    D.hs_merra = ensure_numeric_column(D.hs_merra, height(D));
    D.hs_merra(:) = NaN;
end

if ~ismember("t_merra", string(D.Properties.VariableNames))
    D.t_merra = NaT(height(D),1);
else
    D.t_merra = ensure_datetime_column(D.t_merra, height(D));
    D.t_merra(:) = NaT;
end

valid_match = ~isnat(D.date) & ~isnan(D.lat) & ~isnan(D.lon);
if ~any(valid_match)
    error('No valid rows in D with non-missing date, lat, and lon.');
end

disp('Reading ERA5 time axis...')
time_raw_era5 = ncread(ncfile_era5, 'time');
t_raw_era5 = read_nc_time(ncfile_era5, 'time', time_raw_era5);
t_day_era5 = dateshift(t_raw_era5, 'start', 'day');

disp('Reading grid...')
x = ncread(ncfile_era5, 'x');
y = ncread(ncfile_era5, 'y');

[X,Y] = ndgrid(x,y);   % correct for snod(x,y,time)
[lat_grid, lon_grid] = projinv(projcrs(3408), X, Y);
lon_grid = wrapTo180_local(lon_grid);

nx = numel(x);
ny = numel(y);

disp('Checking ERA5 dimension order...')
dim_order_era5 = get_var_dim_order(ncfile_era5, 'snod');
if ~isequal(dim_order_era5, ["x","y","time"])
    error('ERA5 snod dimension order is %s, expected x,y,time.', strjoin(dim_order_era5, ','));
end

disp('Reading MERRA time axes...')
nfiles_merra = numel(ncfiles_merra);
F = struct([]);

for k = 1:nfiles_merra
    F(k).file = ncfiles_merra{k};
    F(k).time_raw = ncread(F(k).file, 'time');
    F(k).t_raw = read_nc_time(F(k).file, 'time', F(k).time_raw);
    F(k).t_day = dateshift(F(k).t_raw, 'start', 'day');

    dim_order_merra = get_var_dim_order(F(k).file, 'snod');
    if ~isequal(dim_order_merra, ["x","y","time"])
        error('MERRA file %s has snod dimension order %s, expected x,y,time.', ...
            F(k).file, strjoin(dim_order_merra, ','));
    end
end

disp('Computing hs_smlg, t_smlg, hs_merra and t_merra using +/-2 days and 5x5 cells...')
idx_valid = find(valid_match);

for kk = 1:numel(idx_valid)
    i = idx_valid(kk);

    dist_all_km = distance_km_haversine(D.lat(i), D.lon(i), lat_grid, lon_grid);
    [~, ig] = min(dist_all_km(:));
    [ix0, iy0] = ind2sub(size(lat_grid), ig);

    ix = max(1, ix0-grid_half):min(nx, ix0+grid_half);
    iy = max(1, iy0-grid_half):min(ny, iy0+grid_half);

    obs_day = dateshift(D.date(i), 'start', 'day');

    % ---------------- ERA5 ----------------
    vals = [];
    times_used = NaT(0,1);

    if D.date(i) <= cutoff_date_era5
        target_dates = obs_day;
    else
        start_year = year(cutoff_date_era5) - n_clim_years + 1;
        end_year   = year(cutoff_date_era5);
        target_dates = NaT(0,1);
        for yy = start_year:end_year
            td = safe_datetime(yy, month(obs_day), day(obs_day));
            if ~isnat(td)
                target_dates(end+1,1) = td;
            end
        end
    end

    for jj = 1:numel(target_dates)
        t1 = target_dates(jj) - days(time_half);
        t2 = target_dates(jj) + days(time_half);
        ind_t = find(t_day_era5 >= t1 & t_day_era5 <= t2);

        if isempty(ind_t)
            continue
        end

        start = [ix(1), iy(1), ind_t(1)];
        count = [numel(ix), numel(iy), numel(ind_t)];
        sub = ncread(ncfile_era5, 'snod', start, count);
        sub = single_to_nan(sub);

        vals = [vals; sub(:)]; %#ok<AGROW>

        tmp_times = repmat(t_raw_era5(ind_t(:))', numel(ix)*numel(iy), 1);
        times_used = [times_used; tmp_times(:)]; %#ok<AGROW>
    end

    good = ~isnan(vals);
    if any(good)
        vals_good = vals(good);
        times_good = times_used(good);

        D.hs_smlg(i) = mean(vals_good, 'omitnan');

        if D.date(i) <= cutoff_date_era5
            [~, jt] = min(abs(times_good - D.date(i)));
            D.t_smlg(i) = times_good(jt);
        else
            D.t_smlg(i) = max(times_good);
        end
    end

    % ---------------- MERRA ----------------
    vals = [];
    times_used = NaT(0,1);

    if D.date(i) <= cutoff_date_merra
        target_dates = obs_day;
    else
        start_year = year(cutoff_date_merra) - n_clim_years + 1;
        end_year   = year(cutoff_date_merra);
        target_dates = NaT(0,1);
        for yy = start_year:end_year
            td = safe_datetime(yy, month(obs_day), day(obs_day));
            if ~isnat(td)
                target_dates(end+1,1) = td;
            end
        end
    end

    for jj = 1:numel(target_dates)
        t1 = target_dates(jj) - days(time_half);
        t2 = target_dates(jj) + days(time_half);

        for k = 1:nfiles_merra
            ind_t = find(F(k).t_day >= t1 & F(k).t_day <= t2);

            if isempty(ind_t)
                continue
            end

            start = [ix(1), iy(1), ind_t(1)];
            count = [numel(ix), numel(iy), numel(ind_t)];
            sub = ncread(F(k).file, 'snod', start, count);
            sub = single_to_nan(sub);

            vals = [vals; sub(:)]; %#ok<AGROW>

            tmp_times = repmat(F(k).t_raw(ind_t(:))', numel(ix)*numel(iy), 1);
            times_used = [times_used; tmp_times(:)]; %#ok<AGROW>
        end
    end

    good = ~isnan(vals);
    if any(good)
        vals_good = vals(good);
        times_good = times_used(good);

        D.hs_merra(i) = mean(vals_good, 'omitnan');

        if D.date(i) <= cutoff_date_merra
            [~, jt] = min(abs(times_good - D.date(i)));
            D.t_merra(i) = times_good(jt);
        else
            D.t_merra(i) = max(times_good);
        end
    end

    if mod(kk,100) == 0 || kk == numel(idx_valid)
        fprintf('Processed %d / %d rows\n', kk, numel(idx_valid));
    end
end

save(outputFile, 'D')
fprintf('Added ERA5 and MERRA2 snow products.\n')
fprintf('Saved updated dataset to:\n%s\n', ...
    outputFile)

%% Helpers

function x = ensure_numeric_column(x, n)
    if isempty(x)
        x = nan(n,1);
        return
    end
    x = x(:);
    if numel(x) ~= n
        error('Numeric column has inconsistent length.');
    end
    if ~isnumeric(x)
        error('Expected a numeric column.');
    end
end

function t = ensure_datetime_column(t, n)
    if isempty(t)
        t = NaT(n,1);
        return
    end
    t = t(:);
    if isdatetime(t)
        if numel(t) ~= n
            error('Datetime column has inconsistent length.');
        end
        return
    end
    try
        t = datetime(t);
    catch
        error('Could not convert column to datetime.');
    end
    if numel(t) ~= n
        error('Datetime column has inconsistent length.');
    end
end

function lon = wrapTo180_local(lon)
    lon = mod(lon + 180, 360) - 180;
end

function t = safe_datetime(y, m, d)
    try
        t = datetime(y,m,d);
    catch
        t = NaT;
    end
end

function t = read_nc_time(project,varname,time_values)
    units = '';
    try
        units = ncreadatt(project,varname,'units');
    catch
    end

    if isstring(units) || ischar(units)
        units = char(units);
        tok = regexp(units,'(\w+)\s+since\s+(\d{1,4})-(\d{1,2})-(\d{1,2})(?:[ T](\d{1,2}:\d{1,2}:\d{1,2}))?','tokens','once');
        if ~isempty(tok)
            base_unit = lower(tok{1});
            yyyy = str2double(tok{2});
            mm   = str2double(tok{3});
            dd   = str2double(tok{4});

            if numel(tok) >= 5 && ~isempty(tok{5})
                base_time = tok{5};
            else
                base_time = '00:00:00';
            end

            t0 = datetime(sprintf('%04d-%02d-%02d %s',yyyy,mm,dd,base_time), ...
                'InputFormat','yyyy-MM-dd HH:mm:ss');

            switch base_unit
                case {'hour','hours'}
                    t = t0 + hours(time_values);
                case {'day','days'}
                    t = t0 + days(time_values);
                case {'minute','minutes'}
                    t = t0 + minutes(time_values);
                case {'second','seconds'}
                    t = t0 + seconds(time_values);
                otherwise
                    error('Unsupported time unit "%s" in %s.', base_unit, project);
            end
            return
        end
    end

    t0 = datetime(1,1,1,0,0,0);
    t = t0 + hours(time_values);
end

function dim_order = get_var_dim_order(ncfile, varname)
    info = ncinfo(ncfile, varname);
    dim_order = string({info.Dimensions.Name});
end

function A = single_to_nan(A)
    A = double(A);
    A(A <= -9990) = NaN;
end

function d = distance_km_haversine(lat1, lon1, lat2, lon2)
    R = 6371.0088; % mean Earth radius, km

    lat1 = double(lat1);
    lon1 = double(lon1);
    lat2 = double(lat2);
    lon2 = double(lon2);

    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    dlat = lat2 - lat1;
    dlon = mod((lon2 - lon1) + pi, 2*pi) - pi;

    a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
    c = 2 .* atan2(sqrt(a), sqrt(1-a));

    d = R .* c;
end

end