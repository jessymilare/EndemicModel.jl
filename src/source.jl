# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

"""
A data source is supposed to be imported, exported or both.
Subtypes representing importable data source must implement method `import_data`.
Subtypes representing exportable data source must implement method `export_data`.
The method `pathof` should return the path to the database file or directory.
"""
abstract type AbstractDataSource end
abstract type AbstractDataFile <: AbstractDataSource end

const DataSourceDesignator = Union{AbstractDataSource,Symbol}

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
    @debug "Importing data." sources
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
        data =
            import_data(source; cache_data_type = nothing, force_update = true, kwargs...)
        result[sourcesym] = data
    end
    if cache_data_type != nothing
        export_data(cache_data_type, result; kwargs...)
    end
    @debug "Data imported." sources _debuginfo(result)
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
    data = import_data(source; cache_data_type = nothing, force_update = true, kwargs...)
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
    ext::Union{Nothing,AbstractString} = data_extension(option(:cache_data_type));
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

function cache_path(sourcetype::Type{T}; kwargs...) where {T<:AbstractDataFile}
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
_dateformat_aux(per::Minute) =
    Dates.value(per) >= 60 ? dateformat"yyyy-mm-dd HHh" : dateformat"yyyy-mm-dd HHhMM\m"
_dateformat_aux(per::Second) = Dates.value(per) >= 60 ? dateformat"yyyy-mm-dd HHhMM\m" :
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
    @debug "Downloading data." url path
    download(url, path)
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

function import_data(source::CsvPath; simple::Bool = false, kwargs...)
    path = pathof(source)
    @argcheck exists(path)
    if isfile(path)
        @debug "Importing DataFrame from CSV file." path
        data = csv_read(path; copycols = true)
    else
        @debug "Importing Dict from CSV directory." path
        data = DataDict(
            # Parameter `simple` don't need to be passed on
            Symbol(filename(elt)) => import_data(CsvPath(elt); kwargs...)
            for elt ∈ readpath(path)
        )
    end
    if simple
        data = simplify(data; kwargs...)
    end
    data
end

function export_data(
    source::CsvPath,
    data::AbstractDataFrame;
    pretty::Bool = false,
    kwargs...,
)
    destiny = pathof(source)
    @debug "Exporting DataFrame to CSV file." destiny _debuginfo(data)
    if pretty
        data = prettify(data; kwargs...)
    end
    exists(destiny) && rm(destiny; recursive = true)
    csv_write(destiny, data; transform = (col, val) -> something(val, missing))
    destiny
end

function export_data(source::CsvPath, data::AbstractDict; pretty::Bool = false, kwargs...)
    destiny = pathof(source)
    @debug "Exporting Dict to CSV directory." destiny _debuginfo(data)
    if pretty
        data = prettify(data; kwargs...)
    end
    exists(destiny) && rm(destiny; recursive = true, force = true)
    mkdir(destiny; recursive = true)
    for (fname, fdata) ∈ data
        fname = replace(string(fname), '/' => '_')
        fpath = join(destiny, Path(fname * ".csv"))
        export_data(CsvPath(fpath), fdata; pretty = pretty, kwargs...)
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

function import_data(source::OdsPath; simple::Bool = false, kwargs...)
    path = pathof(source)
    @argcheck exists(path)
    if isfile(path)
        data = DataDict(Symbol(k) => v for (k, v) ∈ ods_readall(path))
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
    if simple
        data = simplify(data; kwargs...)
    end
end

function export_data(
    source::OdsPath,
    data::AbstractDataFrame;
    pretty::Bool = false,
    table_name = option(:cache_table_name),
    kwargs...,
)
    if pretty
        data = prettify(data; kwargs...)
    end
    export_data(
        source,
        DataFrameDict(Symbol(table_name) => data);
        pretty = pretty,
        kwargs...,
    )
end

function export_data(
    source::OdsPath,
    data::DataFrameDict;
    pretty::Bool = false,
    kwargs...,
)
    destiny = pathof(source)
    rm(destiny; recursive = true, force = true)
    sheetdict = SortedDict{Tuple,DataFrame}()
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

function export_data(
    source::OdsPath,
    data::AbstractDataDict;
    pretty::Bool = false,
    kwargs...,
)
    if pretty
        data = prettify(data; kwargs...)
    end
    destiny = pathof(source)
    rm(destiny; recursive = true, force = true)
    subdata = filter(p -> p[2] isa AbstractDict, data)
    dataframes = filter(p -> !(p[2] isa AbstractDict), data)
    dataframes = Dict{Symbol,DataFrame}(k => DataFrame(v) for (k, v) ∈ dataframes)
    if isempty(subdata)
        dfsource = source
    else
        mkdir(destiny; recursive = true)
        dfsource = OdsPath(join(destiny, p"_DEFAULT_.ods"))

        for (fname, fdata) ∈ subdata
            fname = replace(string(fname), '/' => '_')
            fpath = join(destiny, Path(fname * ".ods"))
            export_data(OdsPath(fpath), fdata; pretty = pretty, kwargs...)
        end
    end
    export_data(dfsource::OdsPath, dataframes; kwargs...)
    destiny
end
