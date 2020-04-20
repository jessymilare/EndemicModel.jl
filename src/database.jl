abstract type AbstractDatabase{D <: AbstractDict} end
abstract type AbstractSubDatabase{D <: AbstractDict} <: AbstractDatabase{D} end

sourcesof(database::AbstractDatabase) = database.sources
sourcesof!(database::AbstractDatabase, value) = (database.sources = value)
kwargsof(database::AbstractDatabase) = database.kwargs
kwargsof!(database::AbstractDatabase, value) = (database.kwargs = value)
computeof(database::AbstractDatabase) = database.compute
computeof!(database::AbstractDatabase, value) = (database.compute = value)
dataof(database::AbstractDatabase) = database.data
dataof!(database::AbstractDatabase, value) = (database.data = value)

function Base.summary(database::AbstractDatabase{D}) where {D}
    data::D = dataof(database)
    len = length(data)
    nsubdata = count(v -> v isa D, values(data))
    ntables = len - nsubdata
    type = typeof(database).name
    "$(len)-element $(type) [$(ntables) table(s), $(nsubdata) subdatabase(s)]"
end

mutable struct Database{D <: AbstractDict} <: AbstractDatabase{D}
    sources::Vector{Symbol}
    kwargs::Dict{Symbol, Any}
    compute::Vector{Function}
    data::D
    function Database{D}(
        sources::Vector{Symbol},
        kwargs::Dict{Symbol, Any},
        compute::Vector{Function} = Function[],
    ) where {D <: AbstractDict}
        db = new(sources, kwargs, compute, D())
        data = import_data(db; kwargs...)
        eltdata = length(data) == 1 ? collect(values(data))[1] : data
        db.data = eltdata isa AbstractDict ? eltdata : data
        db
    end
end

rootof(database::Database) = database

function Database{D}(
    sources::Vector{Symbol},
    kwargs,
    compute::Vector{Function} = Function[],
) where {D <: AbstractDict}
    Database{D}(sources, Dict{Symbol, Any}(kwargs), compute)
end

Database(args...) = Database{Dict}(args...)

function Database{D}(
    sources::Vector{Symbol},
    compute::Vector{Function} = Function[];
    kwargs...,
) where {D <: AbstractDict}
    Database{D}(sources, Dict(kwargs), compute)
end

function Database{D}(source::Symbol, args...; kwargs...) where {D}
    Database{D}(Vector{Symbol}([source]), args...; kwargs...)
end

function show(io::IO, database::Database{D}) where {D}
    print(io, summary(database), "\n")
end

function show(io::IO, ::MIME"text/plain", database::Database{D}) where {D}
    show(io, database)
end

function show(database::Database{D}) where {D}
    show(stdout, database)
end

mutable struct SubDatabase{D} <: AbstractSubDatabase{D}
    root::Database
    datakeys::Vector{Symbol}
    data::D
    function SubDatabase{D}(
        root::Database,
        datakeys::Vector{Symbol},
        data::D;
        kwargs...,
    ) where {D <: AbstractDict}
        new(root, datakeys, data)
    end
end

sourcesof(database::SubDatabase) = sourcesof(database.root)
sourcesof!(database::SubDatabase, value) = sourcesof!(database.root, value)
kwargsof(database::SubDatabase) = kwargsof(database.root)
kwargsof!(database::SubDatabase, value) = kwargsof!(database.root, value)
computeof(database::SubDatabase) = computeof(database.root)
computeof!(database::SubDatabase, value) = computeof!(database.root, value)
function dataof(database::SubDatabase)
    pdata = dataof(root)
    dataof(database.root)
end
dataof!(database::SubDatabase, value) = dataof!(database.root, value)

rootof(database::SubDatabase) = database.root

function SubDatabase{D}(
    root::Database,
    datakeys::Vector{Symbol};
    kwargs...,
) where {D <: AbstractDict}
    data = dataof(root)
    for k ∈ datakeys
        data = data[k]
    end
    SubDatabase{D}(root, datakeys, data)
end

function SubDatabase{D}(
    parent::SubDatabase{D},
    datakeys::Vector{Symbol};
    kwargs...,
) where {D <: AbstractDict}
    data = dataof(parent)
    for k ∈ datakeys
        data = data[k]
    end
    pdatakeys = datakeysof(parent)
    ndatakeys = pdatakeys ∪ datakeys
    root = rootof(parent)
    SubDatabase{D}(root, ndatakeys, data)
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
    data = dataof(database)
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
    (_, updated) = create_table!(dataof(database), tblname, tablesyms, func; kwargs...)
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
        for (vcols, gfunc) ∈ cols_and_funcs
            append!(combargs, vcols .=> vcols .=> gfunc)
        end
        groupdf = groupby(df, keycols; sort = true)
        groupdict = Dict{Symbol, Any}()
        for (gkey, df) ∈ pairs(groupdf)
            key = Symbol(gkey[1], map(k -> "_" * string(k), Tuple(gkey)[2:end])...)
            key = ensure_unique(groupdict, key)
            groupdict[key] = DataFrame(df)
        end
        combdf = combine(groupdf; combargs...)
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
        dataof(database),
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

function import_data(database::AbstractDatabase; kwargs...)
    imp_kwargs = merge(kwargsof(database), kwargs)
    data = import_data(sourcesof(database); imp_kwargs...)
    funcs = computeof(database)

    while any(func(data) for func ∈ funcs)
    end
    data
end

function import_data!(database::AbstractDatabase; kwargs...)
    dataof!(database, import_data(database; kwargs...))
end

function database_path(
    ext::Union{Nothing, AbstractString} = data_extension(option(:database_data_type));
    database_directory = option(:database_directory),
    database_filename = option(:database_filename),
    kwargs...,
)
    join(Path(database_directory), Path(database_filename) * "." * ext)
end

function database_path(sourcesym::Symbol; kwargs...)
    source = DATA_SOURCES[sourcesym]
    database_path(data_extension(source); kwargs...)
end

function export_data(source::AbstractDataSource, database::AbstractDatabase; kwargs...)
    data = dataof(database)
    export_data(source, data; kwargs...)
end

function export_data(
    database::AbstractDatabase;
    database_data_type = option(:database_data_type),
    database_directory = option(:database_directory),
    database_filename = option(:database_filename),
    kwargs...,
)
    output = database_path(
        database_data_type;
        database_directory = database_directory,
        database_filename = database_filename,
    )
    export_data(database_data_type, database; output = output, kwargs...)
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

@defgroup br_total(per_state) (date,) [(confirmed, deaths) => sum]

@defcolumn closed(recovered, deaths) recovered .+ deaths
@defcolumn active(confirmed, closed) confirmed .- closed
@defgroup world2(world) (country, date) [
    (latitude, longitude) => mean,
    (confirmed, recovered, deaths, closed, active) => sum,
]

@deftable world3(world2, world_population, total_tests) begin
    df = join(world2, world_population; on = :country, kind = :left)
    gdf = join(df, total_tests; on = [:country, :date], kind = :left)
    keys = [:country, :date, :latitude, :longitude]
    vals1 = [:estimated_population, :total_tests, :test_kind]
    vals2 = [:confirmed, :recovered, :deaths, :closed, :active]
    select(gdf, keys ∪ vals1 ∪ vals2)
end

@defgroup per_country(world3) (country,) [
    (date, latitude, longitude) => identity,
    (estimated_population, total_tests, test_kind) => identity,
    (confirmed, recovered, deaths, closed, active) => identity,
]

@defgroup total(per_country) (date,) [
    (latitude, longitude) => mean,
    (estimated_population, total_tests) => sum,
    (test_kind,) => maximum,
    (confirmed, recovered, deaths, closed, active) => sum,
]

@defcolumn μ_closed(deaths, closed) 100.0 .* deaths ./ closed
@defcolumn μ_confirmed(deaths, confirmed) 100.0 .* deaths ./ confirmed

@defcolumn confirmed_per_1mi(confirmed, estimated_population) begin
    1.0e6 .* confirmed ./ estimated_population
end
@defcolumn deaths_per_1mi(deaths, estimated_population) begin
    1.0e6 .* deaths ./ estimated_population
end
@defcolumn recovered_per_1mi(recovered, estimated_population) begin
    1.0e6 .* recovered ./ estimated_population
end
@defcolumn closed_per_1mi(closed, estimated_population) begin
    1.0e6 .* closed ./ estimated_population
end
@defcolumn active_per_1mi(active, estimated_population) begin
    1.0e6 .* active ./ estimated_population
end

@defcolumn tests_per_confirmed(total_tests, confirmed) 1.0 .* total_tests ./ confirmed
@defcolumn tests_per_1mi(total_tests, estimated_population) begin
    1.0e6 .* total_tests ./ estimated_population
end

function copy_tables_to_csse(data::AbstractDict)
    pop = get(data, :world_population, nothing)
    tests = get(data, :total_tests, nothing)
    csse = get(data, :csse, nothing)
    if csse != nothing
        updated = false
        if pop != nothing && !haskey(csse, :world_population)
            csse[:world_population] = pop
            updated = true
        end
        if tests != nothing && !haskey(csse, :total_tests)
            csse[:total_tests] = tests
            true
        else
            updated
        end
    else
        for (_, subdata) ∈ data
            if subdata isa AbstractDict && copy_tables_to_csse(subdata)
                return true
            end
        end
        false
    end
end

function covid_19_database(
    sources = [:csse, :brasil_io, :world_population, :total_tests];
    kwargs...,
)
    funcs = [
        column_closed,
        column_active,
        group_world2,
        group_br_total,
        copy_tables_to_csse,
        table_world3,
        group_per_country,
        group_total,
        column_μ_closed,
        column_μ_confirmed,
        column_confirmed_per_1mi,
        column_deaths_per_1mi,
        column_recovered_per_1mi,
        column_closed_per_1mi,
        column_active_per_1mi,
        column_tests_per_confirmed,
        column_tests_per_1mi,
    ]
    Database{Dict{Symbol, Any}}(sources, kwargs, funcs)
end
