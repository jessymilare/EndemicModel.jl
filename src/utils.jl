# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

const PathDesignator = Union{AbstractString,AbstractPath}

const Float = Float64
const OptFloat = Union{Missing,Float64}
const OptInt = Union{Missing,Int}

const DataDict = Dict{Symbol,Any}
const OrderedDataDict = OrderedDict{Symbol,Any}
const SortedDataDict = SortedDict{Symbol,Any}
const WeakKeyDataDict = WeakKeyDict{Symbol,Any}
const DataFrameDict = Dict{Symbol,DataFrame}
const AbstractDataDict = AbstractDict{Symbol,Any}

datadict(data::AbstractDict) = DataDict(data)
datadict(data::OrderedDict) = OrderedDataDict(data)
datadict(data::SortedDict) = SortedDataDict(data)
datadict(data::WeakKeyDict) = WeakKeyDataDict(data)

function csv_read(file::PathDesignator; kwargs...)::DataFrame
    CSV.read(string(file), DataFrame; kwargs...)
end
function csv_read(io; kwargs...)::DataFrame
    CSV.read(io, DataFrame; kwargs...)
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
function ods_readall(file::PathDesignator; kwargs...)::Dict{String,DataFrame}
    OdsIO.ods_readall(string(file); innerType = "DataFrame", kwargs...)
end
function ods_readall(io; kwargs...)::Dict{String,DataFrame}
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
    ensure_unique(Symbol[sym = obj isa Symbol ? obj : Symbol(string(obj)) for obj ∈ vec])
end

function subdict(dict::AbstractDict, keys)
    if keys(dict) ⊆ keys
        dict
    else
        filter((k, v) -> k ∈ keys, dict)
    end
end

function get_input(
    table::AbstractString;
    directory = option(:parameter_directory),
    kwargs...,
)
    file = joinpath(Path(directory), Path(table * ".csv"))
    csv_read(file; kwargs...)
end

function get_parameters(
    table::AbstractString;
    directory = option(:parameter_directory),
    kwargs...,
)
    df = get_input(table; directory = directory, kwargs...)
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

function prettify(obj::Symbol; case = titlecase, replace_pairs = ["_" => " "], kwargs...)
    if isnothing(case) && isnothing(replace_pairs)
        obj
    else
        Symbol(prettify(string(obj); case = case, replace_pairs = replace_pairs, kwargs...))
    end
end

function prettify(obj::AbstractString; case = nothing, replace_pairs = nothing, kwargs...)
    !isnothing(replace_pairs) && (obj = replace(obj, replace_pairs...))
    if !isnothing(case) && all(c -> Int(c) < 255, obj)
        case(obj)
    else
        obj
    end
end

function prettify(obj::AbstractDataFrame; kwargs...)
    newcols = [prettify(col; kwargs...) for col ∈ eachcol(obj)]
    newnames = prettify(Symbol.(names(obj)); kwargs...)
    DataFrame(newcols, newnames)
end

function prettify(obj::AbstractDict; kwargs...)
    result = empty(obj)
    for (key, value) ∈ obj
        result[prettify(key; kwargs...)] = prettify(value; kwargs...)
    end
    result
end

function prettify(obj::Missing; kwargs...)
    nothing
end

simplify(obj::Any; kwargs...) = obj
function simplify(obj::AbstractVector; kwargs...)
    map(elt -> simplify(elt; kwargs...), obj)
end

function simplify(obj::Symbol; case = lowercase, replace_pairs = [" " => "_"], kwargs...)
    if isnothing(case) && isnothing(replace_pairs)
        obj
    else
        Symbol(simplify(string(obj); case = case, replace_pairs = replace_pairs, kwargs...))
    end
end

function simplify(obj::AbstractString; case = nothing, replace_pairs = nothing, kwargs...)
    !isnothing(replace_pairs) && (obj = replace(obj, replace_pairs...))
    if !isnothing(case) && all(c -> Int(c) < 255, obj)
        case(obj)
    else
        obj
    end
end

function simplify(obj::AbstractDataFrame; kwargs...)
    newcols = [simplify(col; kwargs...) for col ∈ eachcol(obj)]
    newnames = simplify(Symbol.(names(obj)); kwargs...)
    DataFrame(newcols, newnames)
end

function simplify(obj::Nothing; kwargs...)
    missing
end

function simplify(obj::AbstractDict; kwargs...)
    result = empty(obj)
    for (key, value) ∈ obj
        result[simplify(key; kwargs...)] = simplify(value; kwargs...)
    end
    result
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
            jsonval =
                OrderedDict([string(index) => prettify(data[i, index]) for index ∈ values])
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

function _debuginfo(df::AbstractDataFrame)
    summary(df) * " with columns " * "$(Tuple(Symbol.(names(df))))"
end

function _debuginfo(dict::AbstractDict)
    summary(dict) * " with keys " * "$(Tuple(keys(dict)))"
end

_debuginfo(object) = summary(object)

# function maybe_load_language(; force::Bool = false)
#     if Sys.iswindows() && (force || isnothing(get(ENV, "LANGUAGE", nothing)))
#         try
#             txt = read(
#                 `reg query "HKCU\Control Panel\Desktop" /v PreferredUILanguages`,
#                 String,
#             )
#             txt = filter(!isequal('\r'), txt)
#             lines = filter(!isempty, split(txt, '\n'))
#             if lines[1] == "HKEY_CURRENT_USER\\Control Panel\\Desktop"
#                 langs = filter(!isempty, _win_get_lang.(lines[2:end]))
#                 get!(ENV, "LANG", langs[1])
#                 ENV["LANGUAGE"] = join(langs, ':')
#             end
#         catch exc
#             ENV["LANGUAGE"] = "en"
#         end
#     end
#     if force || textdomain() != "endemicmodel"
#         bindtextdomain(
#             "endemicmodel",
#             realpath(joinpath(dirname(Base.pathof(EndemicModel)), "..", "languages")),
#         )
#         textdomain("endemicmodel")
#     end
# end

function _win_get_lang(line)
    tokens = filter(!isempty, split(line, ' '))
    if tokens[1] == "PreferredUILanguages"
        return replace(tokens[end], '-' => '_')
    end
    ""
end
