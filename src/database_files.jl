# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

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

function import_data!(database::AbstractDatabase; kwargs...)
    datadict!(database, import_data(database; kwargs...))
end

function export_data(source::AbstractDataSource, database::AbstractDatabase; kwargs...)
    export_data(source, datadict(database); kwargs...)
end

function export_data(
    database::AbstractDatabase;
    database_data_type = option(:database_data_type),
    database_directory = option(:database_directory),
    database_filename = option(:database_filename),
    model_data_type = option(:model_data_type),
    model_directory = option(:model_directory),
    model_filename = option(:model_filename),
    parameters_data_type = option(:parameters_data_type),
    parameters_subdirectory = option(:parameters_subdirectory),
    parameters_filename = option(:parameters_filename),
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
    parameters_directory = join(model_directory, parameters_subdirectory)
    param_output = database_path(
        parameters_data_type,
        database_directory = parameters_directory,
        database_filename = parameters_filename,
    )

    dbout = export_data(
        database_data_type,
        datadict(database);
        pretty = pretty,
        output = db_output,
        kwargs...,
    )
    modelout = export_data(
        model_data_type,
        modeldata(database);
        pretty = pretty,
        output = model_output,
        kwargs...,
    )
    export_data(
        parameters_data_type,
        paramdict(database);
        pretty = pretty,
        output = param_output,
        kwargs...,
    )
    (dbout, modelout)
end
