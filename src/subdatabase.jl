# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

mutable struct SubDatabase{D<:AbstractDict} <: AbstractDatabase{D}
    parent::AbstractDatabase
    key::Symbol
    datadict::D
    function SubDatabase{D}(
        parent::AbstractDatabase,
        key::Union{Symbol,AbstractString},
        data::D = datadict(parent)[key],
    ) where {D<:AbstractDict}
        new(parent, Symbol(key), data)
    end
end

function SubDatabase{D}(
    parent::SubDatabase{D},
    datakeys::Union{Tuple,AbstractVector},
) where {D<:AbstractDict}
    for key ∈ Symbol.(datakeys)
        data = datadict(parent)
        parent = SubDatabase(parent, key, data[key])
    end
    parent
end

Base.parent(database::SubDatabase) = database.parent

sources(database::SubDatabase) = sources(parent(database))
sources!(database::SubDatabase, value) = sources!(parent(database), value)

default_kwargs(database::SubDatabase) = default_kwargs(parent(database))
default_kwargs!(database::SubDatabase, value) = default_kwargs!(parent(database), value)

computing_functions(database::SubDatabase) = computing_functions(parent(database))
computing_functions!(database::SubDatabase, value) =
    computing_functions!(parent(database), value)

modeldict(database::SubDatabase) = modeldict(parent(database))[database.key]
modeldict!(database::SubDatabase, value) =
    modeldict(parent(database))[database.key] = value

root(database::SubDatabase) = root(parent(database))

function _find_subdata(data::AbstractDict, key)
    result = get(data, key, missing)
    !ismissing(result) && return result

    results = skipmissing(_find_subdata(subdata, key) for subdata ∈ values(data))
    results = collect(results)

    length(results) == 1 && return first(results)
    length(results) > 1 &&
        throw(ErrorException("More than one value found for key $(key)."))
    throw(ErrorException("No value found for key $(key)."))
end

_find_subdata(data::AbstractDatabase, key) = _find_subdata(datadict(data), key)
_find_subdata(data, key) = missing

function Base.getproperty(
    database::AbstractDatabase{D},
    key::Symbol,
) where {D<:AbstractDict}
    if hasfield(typeof(database), key)
        getfield(database, key)
    elseif key == :model
        SubDatabase{D}(database, :model, modeldict(database))
    else
        data = _find_subdata(database, key)
        data = data isa AbstractDict ? SubDatabase{D}(database, key) : data
    end
end

function Base.propertynames(database::AbstractDatabase{D}) where {D<:AbstractDict}
    props = copy(collect(fieldnames(typeof(database))))
    data = datadict(database)
    append!(props, keys(data))
end

function Base.getindex(
    database::AbstractDatabase{D},
    key::Union{Symbol,AbstractString},
    keys::Union{Symbol,AbstractString}...,
) where {D<:AbstractDict}
    database = getproperty(database, Symbol(key))
    for key ∈ keys
        database = getproperty(database, Symbol(key))
    end
    database
end

function Base.setindex!(
    database::AbstractDatabase{D},
    value,
    key::Union{Symbol,AbstractString},
    keys::Union{Symbol,AbstractString}...,
) where {D<:AbstractDict}
    key = Symbol(key)
    if isempty(keys)
        datadict(database)[key] = value
    else
        data, key = get!(datadict(database), key, D()), Symbol(pop!(keys))
        while !isempty(keys)
            data, key = get!(data, key, D()), Symbol(pop!(keys))
        end
        data[key] = value
    end
end

function Base.get(database::AbstractDatabase, key::Union{Symbol,AbstractString}, default)
    get(datadict(database), Symbol(key), default)
end

function Base.get(
    f::Base.Callable,
    database::AbstractDatabase,
    key::Union{Symbol,AbstractString},
)
    get(f, datadict(database), Symbol(key))
end

function Base.get!(database::AbstractDatabase, key::Union{Symbol,AbstractString}, default)
    get!(datadict(database), Symbol(key), default)
end
