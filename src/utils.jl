
const PathDesignator = Union{AbstractString, AbstractPath}

function csv_read(file::PathDesignator; kwargs...)::DataFrame
    CSV.read(string(file); kwargs...)
end
function csv_read(io; kwargs...)::DataFrame
    CSV.read(io; kwargs...)
end
csv_write(file::PathDesignator, table; kwargs...) =
    CSV.write(string(file), table; kwargs...)
csv_write(io, table; kwargs...) = CSV.write(io, table; kwargs...)

function ods_read(file::PathDesignator; kwargs...)::DataFrame
    OdsIO.ods_read(string(file); innerType = "DataFrame", kwargs...)
end
function ods_read(io; kwargs...)::DataFrame
    OdsIO.ods_read(io; innerType = "DataFrame", kwargs...)
end
function ods_readall(file::PathDesignator; kwargs...)::Dict{String, DataFrame}
    OdsIO.ods_readall(string(file); innerType = "DataFrame", kwargs...)
end
function ods_readall(io; kwargs...)::Dict{String, DataFrame}
    OdsIO.ods_readall(io; innerType = "DataFrame", kwargs...)
end
ods_write(file::PathDesignator, data; kwargs...) = OdsIO.ods_write(string(file), data)
ods_write(io, data; kwargs...) = OdsIO.ods_write(io, data)

function ensure_unique(key_collection, key::Symbol, postfix = :_)
    if key ∈ key_collection
        n = 1
        nkey = Symbol(key, postfix, string(n))
        while nkey ∈ key_collection
            n += 1
            nkey = Symbol(key, postfix, string(n))
        end
        nkey
    else
        key
    end
end

function ensure_unique(dict::AbstractDict, key::Symbol, postfix = :_)
    ensure_unique(keys(dict), key, postfix)
end

function ensure_unique(key_collection, ::Nothing, postfix = :_)
    acc = Symbol[key_collection[1]]
    for key::Symbol ∈ key_collection[2:end]
        push!(acc, ensure_unique(acc, key, postfix))
    end
    acc
end

ensure_unique(key_collection) = ensure_unique(key_collection, nothing, :_)

function make_symbols(vec::AbstractVector)
    ensure_unique(Symbol[
        sym = obj isa Symbol ? obj : Symbol(string(obj)) for obj ∈ vec
    ])
end

function subdict(dict::AbstractDict, keys)
    if keys(dict) ⊆ keys
        dict
    else
        filter((k, v) -> k ∈ keys, dict)
    end
end

const _OPTIONS = Dict{Symbol, Any}()

"""
    option(key::Symbol)

Retrieve default option value associated with the given `key`. Returns `default` if no
option is found. See the documentation of `option(::Val{key})` for a given `key`.
"""
function option(key::Symbol, default = missing)
    get(_OPTIONS, key, default)
end

function option(::Val, default = missing)
    default
end

"""
    option!(key::Symbol, value)

Set default option `value` associated with the given `key`. Note that `value` must be
of the type defined with `@defoption`.
"""
function option!(key::Symbol, value)
    option!(Val(key), value)
end

macro defoption(key::Symbol, value, docstring::String, valuetype)
    docstring = """
    option(:$(key))::$(valuetype)

$(docstring)
"""
    key = QuoteNode(key)
    quote
        _OPTIONS[$(key)] = $(value)
        @doc($(docstring), function option(::Val{$(key)})::$(esc(valuetype))
            _OPTIONS[$(key)]
        end)
        function option!(::Val{$(key)}, value::$(esc(valuetype)))
            _OPTIONS[$(key)] = value
        end
    end
end

@defoption(
    cache_directory,
    join(home(), p".cache/EndemicModel/"),
    "Directory where cache files are saved.",
    AbstractPath
)

@defoption(
    parameter_directory,
    join(parent(@__PATH__), p"input/"),
    "Directory to look for parameters.",
    PathDesignator
)

@defoption(
    cache_filename,
    "Sources",
    "Default name of cached data source file.",
    PathDesignator
)

@defoption(cache_data_type, :csv, "Default format of cached data source.", Symbol)

@defoption(
    sources_to_update,
    [:csse, :brasil_io, :world_population, :total_tests],
    "Sources to import from on update.",
    AbstractVector{Symbol}
)

@defoption(cache_nsamples, 7, "Number of cache samples to keep.", Integer)

@defoption(
    cache_table_name,
    "DATABASE",
    "Default name of the table that contains cached imported data.",
    AbstractString
)

@defoption(force_update, false, "If true, force updating data when importing.", Bool)

@defoption(update_period, Day(1), "Frequency of update.", Period)

@defoption(
    database_directory,
    join(parent(@__PATH__), p"database/"),
    "Directory where databases are exported.",
    PathDesignator
)

@defoption(
    database_filename,
    "Database",
    "Default name of exported database file.",
    PathDesignator
)

@defoption(database_data_type, :ods, "Default format of database files.", Symbol)

function get_parameters(
    table::AbstractString;
    parameter_directory = option(:parameter_directory),
)
    file = join(Path(parameter_directory), Path(table) * ".csv")
    df = csv_read(file)
    dfkeys = strip.(df.key)
    dfvalues = strip.(df.value)
    Dict((dfkeys .=> dfvalues)...)
end

function correct_countries!(countries::AbstractVector)
    country_dict = get_parameters("country")
    for (i, country) ∈ enumerate(countries)
        country = strip(country)
        countries[i] = get(country_dict, country, country)
    end
    countries
end

function correct_test_kind!(test_kinds::AbstractVector)
    test_dict = get_parameters("test_kind")
    for (i, test_kind) ∈ enumerate(test_kinds)
        test_kind = strip(test_kind)
        test_kinds[i] = get(test_dict, test_kind, test_kind)
    end
    test_kinds
end

prettify(obj::Any; kwargs...) = obj
function prettify(obj::AbstractVector; kwargs...)
    map(elt -> prettify(elt; kwargs...), obj)
end
prettify(obj::Integer; kwargs...) = obj
function prettify(obj::Real; digits = 2, kwargs...)
    obj < 1.0 && (digits += 1)
    obj < 0.1 && (digits += 1)
    round(obj; digits = digits)
end
prettify(obj::Quantity{<:Integer}; digits = 2, kwargs...) = string(obj)
prettify(obj::Quantity; digits = 2, kwargs...) =
    round(unit(obj), obj; digits = digits + 2)
prettify(obj::Symbol; kwargs...) = Symbol(prettify(string(obj); kwargs...))

function prettify(obj::AbstractString; kwargs...)
    if all(c -> Int(c) < 255, obj)
        titlecase(replace(obj, "_" => " "); strict = false)
    else
        obj
    end
end

function prettify(obj::AbstractDataFrame; kwargs...)
    newcols = [prettify(col; kwargs...) for col ∈ eachcol(obj)]
    newnames = prettify(names(obj); kwargs...)
    DataFrame(newcols, newnames)
end

function prettify(obj::Missing; kwargs...)
    nothing
end

nonmissing(value, default) = ismissing(value) ? default : value

function to_json_dict(
    data::AbstractDataFrame,
    keys = [:state, :city, :date];
    values = [:confirmed, :active, :recovered, :deaths],
    dfvalues = [:estimated_population],
)
    len = nrow(data)
    colnames = names(data)
    key = first(keys)
    restkeys = keys[2:end]
    @argcheck key ∈ colnames

    jsondict = OrderedDict()
    if isempty(restkeys)
        for jsonkey ∈ intersect(dfvalues, colnames)
            jsondict[jsonkey] = data[1, jsonkey]
        end
        values = intersect(values, colnames)
        for i ∈ 1:len
            jsonkey = data[i, key]
            jsonval = OrderedDict(
                [string(index) => prettify(data[i, index]) for index ∈ values],
            )
            jsondict[jsonkey] = jsonval
        end
    else
        for (gkey, gdata) ∈ pairs(groupby(data, key))
            jsondict[gkey[key]] = to_json_dict(gdata, restkeys)
        end
    end
    jsondict
end

function to_json(
    data::AbstractDataFrame,
    keys = [:state, :city, :date];
    values = [:confirmed, :active, :recovered, :deaths],
    dfvalues = [:estimated_population],
)
    JSON.json(to_json_dict(data, keys; values = values, dfvalues = dfvalues), 4)
end
