abstract type AbstractDatabase{D <: AbstractDict} end
abstract type AbstractSubDatabase{D <: AbstractDict} <: AbstractDatabase{D} end

sources(database::AbstractDatabase) = database.sources
sources!(database::AbstractDatabase, value) = database.sources = value

default_kwargs(database::AbstractDatabase) = database.kwargs
default_kwargs!(database::AbstractDatabase, value) = database.kwargs = value

computing_functions(database::AbstractDatabase) = database.computing_functions
computing_functions!(database::AbstractDatabase, value) =
    database.computing_functions = value

datadict(database::AbstractDatabase) = database.datadict
datadict!(database::AbstractDatabase, value) = database.datadict = value

modeldict(database::AbstractDatabase) = database.modeldict
modeldict!(database::AbstractDatabase, value) = database.modeldict = value

paramdict(database::AbstractDatabase) = database.paramdict
paramdict!(database::AbstractDatabase, value) = database.paramdict = value

Base.iterate(database::AbstractDatabase) = iterate(datadict(database))
Base.keys(database::AbstractDatabase) = keys(datadict(database))
Base.values(database::AbstractDatabase) = values(datadict(database))
Base.pairs(database::AbstractDatabase) = pairs(datadict(database))

function Base.summary(database::AbstractDatabase{D}) where {D}
    data::D = datadict(database)
    len = length(data)
    nsubdata = count(v -> v isa D, values(data))
    ntables = len - nsubdata
    type = typeof(database).name
    "$(len)-element $(type) [$(ntables) table(s), $(nsubdata) subdatabase(s)]"
end

mutable struct Database{D <: AbstractDict} <: AbstractDatabase{D}
    sources::Vector{Symbol}
    kwargs::Dict{Symbol, Any}
    computing_functions::Vector{Function}
    datadict::D
    modeldict::D
    paramdict::D
    function Database{D}(
        sources::Vector{Symbol},
        kwargs::Dict{Symbol, Any},
        computing_functions::Vector{Function},
        datadict::Union{D, Nothing},
    ) where {D <: AbstractDict}
        data = something(datadict, D())
        modeldict, paramdict = D(), D()
        db = new(sources, kwargs, computing_functions, data, modeldict, paramdict)
        if isnothing(datadict)
            data = import_data(db; kwargs...)
            eltdata = length(data) == 1 ? collect(values(data))[1] : data
            db.datadict = eltdata isa AbstractDict ? eltdata : data
        end
        db
    end
end

root(database::Database) = database
Base.parent(database::Database) = database

function Database{D}(sources, kwargs, computing_functions) where {D <: AbstractDict}
    sources = collect(Symbol.(sources))
    kwargs, cmp = Dict{Symbol, Any}(kwargs), Vector{Function}(computing_functions)
    Database{D}(sources, kwargs, cmp, nothing)
end

Database(args...) = Database{Dict{Symbol, Any}}(args...)

function Database{D}(
    sources::Union{AbstractVector, Tuple},
    computing_functions::Vector{Function} = Function[];
    kwargs...,
) where {D <: AbstractDict}
    sources = collect(Symbol.(sources))
    Database{D}(sources, Dict{Symbol, Any}(kwargs), computing_functions)
end

function Database{D}(
    source::Union{Symbol, AbstractString},
    args...;
    kwargs...,
) where {D}
    Database{D}(Vector{Symbol}([Symbol(source)]), args...; kwargs...)
end

function Base.show(io::IO, database::AbstractDatabase{D}) where {D}
    print(io, summary(database), "\n")
end

function Base.show(io::IO, ::MIME"text/plain", database::AbstractDatabase{D}) where {D}
    show(io, database)
end

function Base.show(database::AbstractDatabase{D}) where {D}
    show(stdout, database)
end

mutable struct SubDatabase{D <: AbstractDict} <: AbstractDatabase{D}
    parent::AbstractDatabase
    key::Symbol
    datadict::D
    function SubDatabase{D}(
        parent::AbstractDatabase,
        key::Union{Symbol, AbstractString},
        data::D = datadict(parent)[key],
    ) where {D <: AbstractDict}
        new(parent, Symbol(key), data)
    end
end

function SubDatabase{D}(
    parent::SubDatabase{D},
    datakeys::Union{Tuple, AbstractVector},
) where {D <: AbstractDict}
    for key ∈ Symbol.(datakeys)
        data = datadict(parent)
        parent = SubDatabase(parent, key, data[key])
    end
    parent
end

Base.parent(database::SubDatabase) = database.parent

sources(database::SubDatabase) = sources(parent(database))
sources!(database::SubDatabase, value) = sources!(parent(database)value)
default_kwargs(database::SubDatabase) = default_kwargs(parent(database))
default_kwargs!(database::SubDatabase, value) =
    default_kwargs!(parent(database), value)
computing_functions(database::SubDatabase) = computing_functions(parent(database))
computing_functions!(database::SubDatabase, value) =
    computing_functions!(parent(database), value)

root(database::SubDatabase) = root(parent(database))

function Base.getproperty(
    database::AbstractDatabase{D},
    key::Symbol,
) where {D <: AbstractDict}
    if hasfield(typeof(database), key)
        getfield(database, key)
    else
        data = datadict(database)[key]
        data isa AbstractDict ? SubDatabase{D}(database, key) : data
    end
end

function Base.propertynames(database::AbstractDatabase{D}) where {D <: AbstractDict}
    props = copy(collect(fieldnames(typeof(database))))
    append!(props, keys(datadict(database)))
end

function Base.getindex(
    database::AbstractDatabase{D},
    key::Union{Symbol, AbstractString},
    keys::Union{Symbol, AbstractString}...,
) where {D <: AbstractDict}
    database = getproperty(database, Symbol(key))
    for key ∈ keys
        database = getproperty(database, Symbol(key))
    end
    database
end

function Base.setindex!(
    database::AbstractDatabase{D},
    value,
    key::Union{Symbol, AbstractString},
    keys::Union{Symbol, AbstractString}...,
) where {D <: AbstractDict}
    if isempty(keys)
        datadict(database)[key] = value
    else
        data, key = get!(datadict(database), key, D()), pop!(keys)
        while !isempty(keys)
            data, key = get!(data, key, D()), pop!(keys)
        end
        data[key] = value
    end
end

function Base.get(
    database::AbstractDatabase,
    key::Union{Symbol, AbstractString},
    default,
)
    get(datadict(database), Symbol(key), default)
end

function Base.get(
    f::Base.Callable,
    database::AbstractDatabase,
    key::Union{Symbol, AbstractString},
)
    get(f, datadict(database), Symbol(key))
end

function Base.get!(
    database::AbstractDatabase,
    key::Union{Symbol, AbstractString},
    default,
)
    get!(datadict(database), Symbol(key), default)
end

function create_column!(
    df::AbstractDataFrame,
    colname::Symbol,
    columnsyms,
    func::Function;
    kwargs...,
)
    dfcols = names(df)
    if hasproperty(df, colname) || any(col -> !hasproperty(df, col), columnsyms)
        (df, false)
    else
        ind = ncol(df) + 1
        columns = [df[!, col] for col ∈ columnsyms]
        newcol = func(columns...)
        df = insertcols!(df, ind, colname => newcol)
        (df, true)
    end
end

function create_column!(
    data::AbstractDict,
    colname::Symbol,
    columnsyms,
    func::Function;
    kwargs...,
)
    updated = false
    for (_, elt) ∈ data
        (_, elt_updt) = create_column!(elt, colname, columnsyms, func; kwargs...)
        updated = updated || elt_updt
    end
    (data, updated)
end

function create_column!(
    database::AbstractDatabase,
    colname::Symbol,
    columnsyms,
    func::Function;
    kwargs...,
)
    data = datadict(database)
    (_, updated) = create_column!(data, colname, columnsyms, func; kwargs...)
    (database, updated)
end

function create_table!(
    data::AbstractDict,
    tblname::Symbol,
    tablesyms,
    func::Function;
    kwargs...,
)
    updated = false
    for (key, elt) ∈ data
        if elt isa AbstractDict
            (df, elt_updt) = create_table!(elt, tblname, tablesyms, func; kwargs...)
            updated = updated || elt_updt
        end
    end
    if !haskey(data, tblname) && all(tbl -> haskey(data, tbl), tablesyms)
        updated = true
        tables = [data[tbl] for tbl ∈ tablesyms]
        data[tblname] = func(tables...)
    else
        updated = false
    end
    (data, updated)
end

function create_table!(
    database::AbstractDatabase,
    tblname::Symbol,
    tablesyms,
    func::Function;
    kwargs...,
)
    (_, updated) =
        create_table!(datadict(database), tblname, tablesyms, func; kwargs...)
    (database, updated)
end

function group_table!(
    data::AbstractDict,
    tblname::Symbol,
    sourcetbl::Symbol,
    keycols,
    cols_and_funcs::Pair...;
    kwargs...,
)
    if !haskey(data, tblname) && haskey(data, sourcetbl)
        df = data[sourcetbl]
        keycols = Vector{Symbol}(keycols)
        combargs = []
        for (vcol, gfunc) ∈ cols_and_funcs
            append!(combargs, vcol .=> (col -> (; vcol => gfunc(col))))
        end
        groupdf = groupby(df, keycols; sort = true)
        groupdict = Dict{Symbol, Any}()
        for (gkey, df) ∈ pairs(groupdf)
            key = Symbol(gkey[1], map(k -> "_" * string(k), Tuple(gkey)[2:end])...)
            key = ensure_unique(groupdict, key)
            groupdict[key] = DataFrame(df)
        end
        combdf = combine(groupdf) do subdf
            colpairs = []
            for (vcols, gfunc) ∈ cols_and_funcs, vcol ∈ vcols
                push!(colpairs, vcol => gfunc(subdf[!, vcol]))
            end
            DataFrame(; colpairs...)
        end
        tblname_groups = ensure_unique(data, Symbol(tblname, :_group))
        data[tblname_groups] = groupdict
        data[tblname] = combdf
        updated = true
    else
        updated = false
    end
    for (_, elt) ∈ data
        if elt isa AbstractDict
            (df, elt_updt) = group_table!(
                elt,
                tblname,
                sourcetbl,
                keycols,
                cols_and_funcs...;
                kwargs...,
            )
            updated = updated || elt_updt
        end
    end
    (data, updated)
end

function group_table!(
    database::AbstractDatabase,
    tblname::Symbol,
    sourcetbl::Symbol,
    keycols,
    cols_and_funcs::Pair...;
    kwargs...,
)
    (_, updated) = group_table!(
        datadict(database),
        tblname,
        sourcetbl,
        keycols,
        cols_and_funcs...;
        kwargs...,
    )
    (database, updated)
end

function group_table!(data, tblname, sourcetbls, keycols, cols_and_funcs...; kwargs...)
    updated = false
    for sourcetbl ∈ sourcetbls
        (_, updt) = group_table!(
            data,
            tblname,
            sourcetbl,
            keycols,
            cols_and_funcs...;
            kwargs...,
        )
        updated = updated || updt
    end
    (data, updated)
end

function import_data(database::AbstractDatabase; kw...)
    imp_kwargs = merge(default_kwargs(database), kw)
    dt = import_data(sources(database); imp_kwargs...)
    funcs = computing_functions(database)

    while any(func(dt) for func ∈ funcs)
    end
    dt
end

macro defcolumn(fcall::Expr, body::Expr)
    @argcheck fcall.head == :call
    colname = fcall.args[1]
    funname = Symbol(:column_, colname)
    columnnames = fcall.args[2:end]
    args = Tuple(columnnames)
    colsym = QuoteNode(colname)
    quote
        function $(esc(funname))(data)
            func = ($(esc.(args)...),) -> $(esc(body))
            create_column!(data, $(colsym), $(args), func)[2]
        end
    end
end

macro deftable(fcall::Expr, body::Expr)
    @argcheck fcall.head == :call
    tblname = fcall.args[1]
    funname = Symbol(:table_, tblname)
    tablenames = fcall.args[2:end]
    args = Tuple(tablenames)
    tblsym = QuoteNode(tblname)
    quote
        function $(esc(funname))(data)
            func = ($(esc.(args)...),) -> $(esc(body))
            create_table!(data, $(tblsym), $(args), func)[2]
        end
    end
end

macro defgroup(fcall::Expr, keynames::Expr, cols_and_funcs::Expr)
    @argcheck fcall.head == :call
    @argcheck keynames.head ∈ (:tuple, :vect)
    @argcheck cols_and_funcs.head ∈ (:tuple, :vect, :block)

    tblname = fcall.args[1]
    funname = Symbol(:group_, tblname)
    tablesyms = fcall.args[2:end]
    args = Tuple(tablesyms)
    tblsym = QuoteNode(tblname)
    keysyms = [k isa QuoteNode ? k.value : k for k ∈ keynames.args]

    cols_funcs = []
    for pair ∈ cols_and_funcs.args
        if !(pair isa LineNumberNode)
            @argcheck pair isa Expr && pair.head == :call && pair.args[1] == :(=>)
            (cols, func) = filter(elt -> !(elt isa LineNumberNode), pair.args[2:end])
            @argcheck cols isa Expr && cols.head ∈ (:tuple, :vect)

            push!(cols_funcs, :($(Tuple(cols.args)) => $(esc(func))))
        end
    end

    quote
        function $(esc(funname))(data)
            group_table!(data, $(tblsym), $(args), $keysyms, $(cols_funcs...))[2]
        end
    end
end
