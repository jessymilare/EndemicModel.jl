# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

function database_path(
    ext::Union{Nothing,AbstractString} = data_extension(option(:database_data_type));
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

function import_data(
    database::AbstractDatabase;
    database_data_type = option(:database_data_type),
    database_directory = option(:database_directory),
    database_filename = name(database),
    model_data_type = option(:model_data_type),
    model_directory = option(:model_directory),
    model_filename = name(database),
    simple::Bool = true,
    force_update::Bool = false,
    kwargs...,
)
    @debug "Creating database."
    imported = false
    if !force_update
        db_input = database_path(
            database_data_type;
            database_directory = database_directory,
            database_filename = database_filename,
        )
        model_input = database_path(
            model_data_type;
            database_directory = model_directory,
            database_filename = database_filename,
        )
        if (
            exists(db_input) &&
            exists(model_input) &&
            modified(db_input) >= today() &&
            modified(model_input) >= today()
        )
            @debug "Importing data from cache." db_input model_input
            data = import_data(
                database_data_type;
                path = db_input,
                simple = simple,
                kwargs...,
            )
            modeldata = import_data(
                database_data_type;
                path = model_input,
                simple = simple,
                kwargs...,
            )
            model = SEIRModel(data; modeldata = modeldata, kwargs...)
            imported = true
            @debug "Data imported from cache."
        end
    end
    if !imported
        @debug "Importing data from sources."
        imp_kwargs = merge(default_kwargs(database), kwargs...)
        data = import_data(sources(database); imp_kwargs...)
        model = nothing
        funcs = computing_functions(database)

        @debug "Data imported from sources." _debuginfo(data)
        @debug "Computing database."
        while any(func(data) for func ∈ funcs)
        end
    end
    @debug "Database created."
    (data, model)
end

function import_data!(database::AbstractDatabase; kwargs...)
    (data, model) = import_data(database; kwargs...)
    datadict!(database, data)
    !isnothing(model) && modeldict!(database, model)
    (data, model)
end

function export_data(source::AbstractDataSource, database::AbstractDatabase; kwargs...)
    export_data(source, datadict(database); kwargs...)
end

function export_data(
    database::AbstractDatabase;
    database_data_type = option(:database_data_type),
    database_directory = option(:database_directory),
    database_filename = name(database),
    model_data_type = option(:model_data_type),
    model_directory = option(:model_directory),
    model_filename = name(database),
    pretty::Bool = true,
    kwargs...,
)
    db_output = database_path(
        database_data_type;
        database_directory = database_directory,
        database_filename = database_filename,
    )
    model_output = database_path(
        model_data_type;
        database_directory = model_directory,
        database_filename = database_filename,
    )
    data = datadict(database)
    @info "Exporting database." path = db_output _debuginfo(data)
    dbout = export_data(
        database_data_type,
        data;
        pretty = pretty,
        output = db_output,
        kwargs...,
    )
    data = modeldict(database)
    @info "Exporting database models." path = model_output _debuginfo(data)
    modelout = export_data(
        model_data_type,
        data;
        pretty = pretty,
        output = model_output,
        kwargs...,
    )
    (dbout, modelout)
end
