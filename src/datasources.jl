"""
A data source is supposed to be imported, exported or both.
Subtypes representing importable data source must implement method `import_data`.
Subtypes representing exportable data source must implement method `export_data`.
The method `pathof` should return the path to the database file or directory.
"""
abstract type AbstractDataSource end
abstract type AbstractDataFile <: AbstractDataSource end

const DataSourceDesignator = Union{AbstractDataSource, Symbol}
const DataDict = Dict{Symbol, Any}
const DataFrameDict = Dict{Symbol, DataFrame}
const AbstractDataDict = AbstractDict{Symbol, Any}

"""
    pathof(source::AbstractDataFile)

Return the path to the file or directory of `source`. Defaults to `source.path`.
"""
function pathof(source::AbstractDataFile)
    hasproperty(source, :path) ? source.path : nothing
end

"""
    import_data(source; kwargs...)

Import data from `source` and converts it either into a table or a named collection
thereof. `source` can be a `Symbol` which will be searched in `DATA_SOURCES`.
"""
function import_data end

"""
    import_data(sources::AbstractVector, data; kwargs...)

Import each source in `sources` and return a `Dict`:
```
Dict(source => import_data(source, data; kwargs...) for source ∈ sources)
```
"""
function import_data(
    sources::AbstractVector = option(:sources_to_update);
    cache_data_type = option(:cache_data_type),
    force_update::Bool = option(:force_update),
    kwargs...,
)
    @argcheck !isempty(sources)
    source_keys = make_symbols(sources)
    try
        if !force_update && cache_data_type != nothing
            data = import_data(
                cache_data_type;
                cache_data_type = nothing,
                force_update = true,
                kwargs...,
            )
            if !(data isa AbstractDict)
                return data
            elseif source_keys ⊆ keys(data)
                return subdict(data, source_keys)
            end
        end
    catch
    end

    result = DataDict()
    for (source, sourcesym) ∈ zip(sources, source_keys)
        data = import_data(
            source;
            cache_data_type = nothing,
            force_update = true,
            kwargs...,
        )
        result[sourcesym] = data
    end
    if cache_data_type != nothing
        export_data(cache_data_type, result; kwargs...)
    end
    result
end

function import_data(
    sourcesym::Symbol;
    path = nothing,
    cache_data_type = option(:cache_data_type),
    force_update::Bool = option(:force_update),
    kwargs...,
)
    try
        if !force_update && cache_data_type != nothing
            data = import_data(
                cache_data_type;
                cache_data_type = nothing,
                force_update = false,
                kwargs...,
            )
            if !(data isa AbstractDict)
                return data
            elseif haskey(data, sourcesym)
                return data[sourcesym]
            end
        end
    catch
    end

    source = DATA_SOURCES[sourcesym](nothing; path = path, kwargs...)
    data =
        import_data(source; cache_data_type = nothing, force_update = true, kwargs...)
    if cache_data_type != nothing
        export_data(cache_data_type, Dict(sourcesym => data); kwargs...)
    end
    data
end

"""
    export_data(source, data; kwargs...)

Export `data` to `source`; `data` must be a table or a named collection thereof.
If `source` can be a `Symbol`, search it in `DATA_SOURCES`.
"""
function export_data end

"""
    export_data(sources::AbstractVector, data; kwargs...)

Export each source in `sources` and return an `Dict`:
```
Dict(source => export_data(source, data; kwargs...) for source ∈ sources)
```
"""
function export_data(
    sources::AbstractVector,
    data;
    source_keys::AbstractVector{Symbol} = make_symbols(sources),
    kwargs...,
)
    Dict(
        sourcesym => export_data(source; kwargs...)
        for (source, sourcesym) ∈ zip(sources, source_keys)
    )
end

function export_data(sourcesym::Symbol, data; output = nothing, kwargs...)
    source = DATA_SOURCES[sourcesym](nothing; path = output, kwargs...)
    export_data(source, data; kwargs...)
end

function export_data(data; cache_data_type = option(:cache_data_type), kwargs...)
    export_data(cache_data_type, data; kwargs...)
end

"""
    cache_path(
        ext::Union{Nothing, AbstractString} = databaseExtension;
        directory = cache_directory,
        sub = cache_filename,
        nsamples::Integer = cache_nsamples,
        kwargs...,
    )

Returns the default cache file path with the given `ext`ension, or no extension if `ext`
is `nothing`.
"""
function cache_path(
    ext::Union{Nothing, AbstractString} = data_extension(option(:cache_data_type));
    directory = option(:cache_directory),
    update_period::Period = option(:update_period),
    subpath = option(:cache_filename),
    nsamples::Integer = option(:cache_nsamples),
    kwargs...,
)
    @argcheck nsamples >= 1
    (period_directory, lastdir) = _cache_subpath(Path(directory), update_period)
    subpath = ext == nothing ? subpath : subpath * "." * ext
    result = join(period_directory, lastdir, Path(subpath))
    curdir = dirname(result)
    if !exists(curdir)
        mkdir(curdir, recursive = true)
        gc_cache_directory(period_directory, nsamples)
    end
    result
end

function _cache_subpath(directory::AbstractPath, update_period::Period)
    try
        update_period = update_period < Second(1) ? Second(1) : update_period
    catch
    end
    n = Dates.value(update_period)
    timeformat = _dateformat_aux(update_period)
    period_directory = join(directory, Path("Update every " * string(update_period)))
    (perfloor, perceil) = Dates.floorceil(now(), update_period)
    if n == 1
        lastdir = Path(Dates.format(perfloor, timeformat))
    else
        lastdir = Path(
            Dates.format(perfloor, timeformat) *
            " to " *
            Dates.format(perceil, timeformat),
        )
    end
    (period_directory, lastdir)
end

function cache_path(sourcetype::Type{T}; kwargs...) where {T <: AbstractDataFile}
    cache_path(data_extension(T); kwargs...)
end
function cache_path(sourcesym::Symbol; kwargs...)
    source = DATA_SOURCES[sourcesym]
    cache_path(data_extension(source); kwargs...)
end

_dateformat_aux(::Year) = dateformat"yyyy"
_dateformat_aux(per::Month) =
    Dates.value(per) >= 12 ? dateformat"yyyy" : dateformat"yyyy-mm"
_dateformat_aux(per::Week) =
    Dates.value(per) >= 5 ? dateformat"yyyy-mm" : dateformat"yyyy-mm-dd"
_dateformat_aux(per::Day) =
    Dates.value(per) >= 31 ? dateformat"yyyy-mm" : dateformat"yyyy-mm-dd"
_dateformat_aux(per::Hour) =
    Dates.value(per) >= 24 ? dateformat"yyyy-mm-dd" : dateformat"yyyy-mm-dd HHh"
_dateformat_aux(per::Minute) = Dates.value(per) >= 60 ? dateformat"yyyy-mm-dd HHh" :
    dateformat"yyyy-mm-dd HHhMM\m"
_dateformat_aux(per::Second) =
    Dates.value(per) >= 60 ? dateformat"yyyy-mm-dd HHhMM\m" :
    dateformat"yyyy-mm-dd HHhMM\mSS\s"

_dateformat_aux(period::Period) = dateformat"yyyy-mm-dd HHhMM\mSS\s"

"""
    gc_cache_directory(directory::AbstractPath, nsamples::Integer)

Garbage collect `directory`. Each subdirectory in cache `directory` is named by
date and/or time it was saved. This function keeps only the oldest `nsamples`
subdirectories.
"""
function gc_cache_directory(directory::AbstractPath, nsamples::Integer)
    @argcheck nsamples >= 1
    dirlist = readpath(directory)
    ndel = length(dirlist) - nsamples
    if ndel > 0
        todelete = sort!(dirlist; by = basename)[1:ndel]
        foreach(d -> rm(d; recursive = true), todelete)
    end
end

tableKeys = [:country, :state, :date]
tableColumns =
    [:country, :state, :latitude, :longitude, :date, :confirmed, :recovered, :deaths]

"""
Auxiliar structure to handle downloads.
"""
struct DownloadSource <: AbstractDataSource
    url::AbstractString
    function DownloadSource(url; kwargs...)
        new(string(url))
    end
end

function import_data(data::DownloadSource; kwargs...)
    url = data.url
    ext = split(url, '.')[end]
    ext = all(isletter, ext) ? ext : nothing
    subpath = p"_DOWNLOADS_/" / string(hash(url))
    path = cache_path(ext; subpath = subpath, kwargs...)
    download(url, path)
end

# CSSE COVID-19 Repository

const CSSE_COLUMNS = (;
    Symbol("Country/Region") => :country,
    Symbol("Province/State") => :state,
    Symbol("Lat") => :latitude,
    Symbol("Long") => :longitude,
)

function _csse_url(tablename::String)
    "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/" *
    "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_" *
    tablename *
    "_global.csv"
end

const CSSE_TABLES = (
    confirmed = _csse_url("confirmed"),
    deaths = _csse_url("deaths"),
    recovered = _csse_url("recovered"),
)

"""
Novel Coronavirus (COVID-19) Cases, provided by JHU CSSE
https://github.com/CSSEGISandData/COVID-19
https://systems.jhu.edu/research/public-health/ncov/
"""
struct CsseSource <: AbstractDataSource
    function CsseSource()
        new()
    end
end
CsseSource(::Nothing; kwargs...) = CsseSource()

DATA_SOURCES[:csse] = CsseSource

function _import_data(source::CsseSource; kwargs...)
    result = DataDict()
    for (dfsym, url) ∈ pairs(CSSE_TABLES)
        file = import_data(DownloadSource(url))
        df = csv_read(file; copycols = true)

        newcolnames = [
            if all(c -> c ∈ "0123456789/", string(colname))
                date_value = Date(string(colname), "m/d/y") + Year(2000)
                Symbol(string(date_value))
            else
                CSSE_COLUMNS[colname]
            end for colname ∈ names(df)
        ]
        rename!(df, newcolnames)
        df[!, :country] = correct_countries!(df.country)
        result[dfsym] = df
    end
    result
end

function import_data(source::CsseSource; kwargs...)
    raw_data = _import_data(source)
    key_columns = [:country, :state, :latitude, :longitude]
    dataframes = []

    for (main_col_name, df) ∈ raw_data
        dates = []
        date_syms = []
        for colsym ∈ names(df)
            colstr = string(colsym)
            if all(c -> c ∈ "0123456789-", colstr)
                push!(dates, Date(colstr))
                push!(date_syms, colsym)
            end
        end
        datadf = by(df, key_columns) do subdf
            col_values = [sum(subdf[!, date_sym]) for date_sym ∈ date_syms]
            DataFrame(; date = dates, main_col_name => col_values)
        end
        push!(dataframes, datadf)
    end
    parsedata = join(dataframes...; on = key_columns ∪ [:date]) |> DataFrame
    sort!(parsedata, [:country, :state, :date])
    DataDict(:world => parsedata)
end

# World Bank Population database

const WORLD_POP_WB_URL = "http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv"
const WORLD_POP_YEAR_0 = 2000

"""
Total population of each country, extracted from https://data.worldbank.org/
"""
struct WorldPopulationWB <: AbstractDataSource
    function WorldPopulationWB()
        new()
    end
end
WorldPopulationWB(::Nothing; kwargs...) = WorldPopulationWB()

DATA_SOURCES[:world_population] = WorldPopulationWB

function import_data(source::WorldPopulationWB; url = WORLD_POP_WB_URL, kwargs...)
    curyear = year(today())
    years = [string(y) for y ∈ WORLD_POP_YEAR_0:curyear]
    keycols = ["Country Name", "Country Code"]
    col_select = keycols ∪ years

    file = import_data(DownloadSource(url))
    global r = ZipFile.Reader(file)
    zipfiles = r.files
    file = zipfiles[findfirst(f -> f.name[1:10] == "API_SP.POP", zipfiles)]

    data = csv_read(file; copycols = true, header = 5, select = col_select)
    close(file)

    rename!(data, "Country Name" => :country, "Country Code" => :country_code)
    data[!, :country] = correct_countries!(data.country)
    sort!(data, :country)

    ncols = ncol(data)
    for nr ∈ ncols:-1:1
        if !all(data[!, nr] .=== missing)
            if nr < ncols
                select!(data, Not((nr + 1):ncols))
            end
            break
        end
    end
    years = intersect(map(Symbol, years), names(data))
    fyear = years[1]
    lyear = years[end]

    deleterows!(data, data[!, fyear] .=== missing)
    myear = years[findlast(y -> all(data[!, y] .!== missing), years)]

    fyear_int = parse(Int, string(fyear))
    lyear_int = parse(Int, string(lyear))
    myear_int = parse(Int, string(myear))

    n = (lyear_int - fyear_int)
    avg_pop_factor = (data[!, lyear] ./ data[!, fyear]) .^ (1 / n)
    navg_factor = Float64[
        factor === missing ? (data[i, myear] / data[i, fyear])^(1 / n) : factor
        for (i, factor) ∈ enumerate(avg_pop_factor)
    ]
    insertcols!(data, ncol(data) + 1; avg_pop_factor = navg_factor)

    n = (curyear - lyear_int)
    pop_est = (data[!, lyear] .* (data.avg_pop_factor .^ n))
    npop_est = [
        val === missing ? data[i, myear] .* (data[i, :avg_pop_factor] .^ n) : val
        for (i, val) ∈ enumerate(pop_est)
    ]
    npop_est = npop_est .|> round .|> Int
    insertcols!(data, ncol(data) + 1; estimated_population = npop_est)

    select(data, :country, :country_code, :avg_pop_factor, :estimated_population)
end

"""
Number of tests per country from Our World in Data.
https://ourworldindata.org/covid-testing
"""
struct TotalTestsOWID <: AbstractDataSource
    function TotalTestsOWID()
        new()
    end
end
TotalTestsOWID(::Nothing; kwargs...) = TotalTestsOWID()

DATA_SOURCES[:total_tests] = TotalTestsOWID

const TESTING_OWID_URL = "https://github.com/owid/covid-19-data/raw/master/public/data/testing/covid-testing-all-observations.csv"

function import_data(source::TotalTestsOWID; url = TESTING_OWID_URL, kwargs...)
    file = import_data(DownloadSource(url))
    col_select = [:Entity, :Date, Symbol("Cumulative total")]
    data::DataFrame = csv_read(file; copycols = true, select = col_select)
    rename!(data, [:country, :date, :total_tests])
    country_split = map(name -> split(name, " - "), data.country)
    data[!, :country] = map(row -> row[1], country_split) |> correct_countries!
    test_kinds = map(row -> row[2], country_split) |> correct_test_kind!
    insertcols!(data, 3, :test_kind => test_kinds)

    rows = collect(eachrow(data))
    nextrows = rows[2:end]
    lastrow = rows[end]
    nlastrow = DataFrame(;
        country = [lastrow.country],
        date = [today()],
        test_kind = [lastrow.test_kind],
        total_tests = [lastrow.total_tests],
    )
    push!(nextrows, nlastrow[1, :])
    for (row, nextrow) ∈ zip(rows, nextrows)
        if row.country == nextrow.country
            dtend = nextrow.date - Day(1)
        else
            dtend = today() - Day(1)
        end
        for dt ∈ (row.date + Day(1)):Day(1):dtend
            append!(
                data,
                DataFrame(;
                    country = [row.country],
                    date = [dt],
                    test_kind = [row.test_kind],
                    total_tests = [row.total_tests],
                ),
            )
        end
    end
    sort!(data, [:country, :date])
    data
end

# Brasil.io covid-19 database

const BRASIL_IO_URL = "https://brasil.io/dataset/covid19/caso?format=csv"

struct BrasilIo <: AbstractDataSource
    function BrasilIo()
        new()
    end
end
BrasilIo(::Nothing; kwargs...) = BrasilIo()

DATA_SOURCES[:brasil_io] = BrasilIo

function import_data(source::BrasilIo; url = BRASIL_IO_URL, kwargs...)
    col_drop = [:confirmed_per_100k_inhabitants, :death_rate, :is_last]
    file = import_data(DownloadSource(url))
    raw_data::DataFrame = csv_read(
        file;
        copycols = true,
        truestrings = ["True"],
        falsestrings = ["False"],
        drop = col_drop,
    )
    raw_data[(raw_data.deaths .=== missing), :deaths] .= 0
    sort!(raw_data, [:state, :city, :date])

    rename!(raw_data, "estimated_population_2019" => :estimated_population)
    per_state = @where(raw_data, :place_type .== "state")
    per_city = @where(raw_data, :place_type .== "city")

    select!(per_state, Not([:city, :place_type]))
    select!(per_city, Not([:place_type]))

    DataDict(:per_state => per_state, :per_city => per_city)
end

# CSV file and directory of files
"""
Path to a CSV file or directory
"""
struct CsvPath <: AbstractDataFile
    path::AbstractPath
end
function CsvPath(::Nothing; path = nothing, kwargs...)
    if path == nothing
        path = cache_path("csv"; kwargs...)
    end
    CsvPath(Path(path))
end

DATA_SOURCES[:csv] = CsvPath
data_extension(::Type{CsvPath}) = "csv"

function import_data(source::CsvPath; kwargs...)
    path = pathof(source)
    @argcheck exists(path)
    if isfile(path)
        csv_read(path; copycols = true)
    else
        DataDict(
            Symbol(filename(elt)) => import_data(CsvPath(elt); kwargs...)
            for elt ∈ readpath(path)
        )
    end
end

function export_data(
    source::CsvPath,
    data::AbstractDataFrame;
    pretty::Bool = false,
    kwargs...,
)
    destiny = pathof(source)
    if pretty
        data = prettify(data)
    end
    exists(destiny) && rm(destiny; recursive = true)
    csv_write(destiny, data; transform = (col, val) -> something(val, missing))
    destiny
end

function export_data(source::CsvPath, data::AbstractDict; kwargs...)
    destiny = pathof(source)
    exists(destiny) && rm(destiny; recursive = true)
    mkdir(destiny; recursive = true)
    for (fname, fdata) ∈ data
        fpath = join(destiny, Path(string(fname) * ".csv"))
        export_data(CsvPath(fpath), fdata; kwargs...)
    end
    destiny
end

struct OdsPath <: AbstractDataFile
    path::AbstractPath
end
function OdsPath(::Nothing; path = nothing, kwargs...)
    if path == nothing
        path = cache_path("ods"; kwargs...)
    end
    OdsPath(Path(path))
end

DATA_SOURCES[:ods] = OdsPath
data_extension(::Type{OdsPath}) = "ods"

function import_data(source::OdsPath; kwargs...)
    path = pathof(source)
    @argcheck exists(path)
    if isfile(path)
        data = ods_readall(path)
        DataDict(Symbol(k) => v for (k, v) ∈ data)
    else
        data = DataDict(
            Symbol(filename(elt)) => import_data(OdsPath(elt); kwargs...)
            for elt ∈ readpath(path)
        )
        dfdict = get(data, :_DEFAULT_, nothing)
        if dfdict != nothing
            delete!(dfdict::Dict, :_DEFAULT_)
            merge!(data, dfdict)
        end
        data
    end
end

function export_data(
    source::OdsPath,
    data::AbstractDataFrame;
    table_name = option(:cache_table_name),
    kwargs...,
)
    export_data(source, DataFrameDict(Symbol(table_name) => data); kwargs...)
end

function export_data(
    source::OdsPath,
    data::DataFrameDict;
    pretty::Bool = false,
    kwargs...,
)
    destiny = pathof(source)
    rm(destiny; recursive = true, force = true)
    sheetdict = SortedDict{Tuple, DataFrame}()
    for (tname, tbl) ∈ data
        tname = string(tname)
        if pretty
            tname = prettify(tname; kwargs...)
            tbl = prettify(tbl; kwargs...)
        end
        push!(sheetdict, (tname, 1, 1) => tbl)
    end
    ods_write(destiny, sheetdict)
    destiny
end

function export_data(source::OdsPath, data::AbstractDataDict; kwargs...)
    destiny = pathof(source)
    rm(destiny; recursive = true, force = true)
    subdata = filter(p -> p[2] isa AbstractDict, data)
    dataframes = filter(p -> !(p[2] isa AbstractDict), data)
    dataframes = Dict{Symbol, DataFrame}(k => DataFrame(v) for (k, v) ∈ dataframes)
    if isempty(subdata)
        dfsource = source
    else
        mkdir(destiny; recursive = true)
        dfsource = OdsPath(join(destiny, p"_DEFAULT_.ods"))

        for (fname, fdata) ∈ subdata
            fpath = join(destiny, Path(string(fname) * ".ods"))
            export_data(OdsPath(fpath), fdata; kwargs...)
        end
    end
    export_data(dfsource::OdsPath, dataframes; kwargs...)
    destiny
end
