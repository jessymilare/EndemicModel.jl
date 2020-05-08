
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

@defoption(database_data_type, :csv, "Default format of database files.", Symbol)

@defoption(
    minimum_confirmed_factor,
    1e-5,
    "Factor for minimum number of confirmed people to be considered in estimates.",
    Real
)
@defoption(
    minimum_plot_factor,
    1e-3,
    "Factor for minimum number of confirmed people to be considered in plots.",
    Real
)