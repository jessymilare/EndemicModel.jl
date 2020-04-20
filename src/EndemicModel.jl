module EndemicModel

import Base.Iterators
import CSV
import OdsIO
import ZipFile

using ArgCheck
using Base
using Base.Iterators
using Dates
using DataFrames
using DataFramesMeta
using DataStructures
using DifferentialEquations
using FilePaths
using Interpolations
using Optim
using Plots
using Statistics
using Unitful

export AbstractDataSource,
    AbstractDataPath,
    CsseSource,
    CsvDirectory,
    CsvPath,
    OdsPath,

    option,
    option!,
    get_parameters,

    DATA_SOURCES,
    import_data,
    export_data,

    Database,

    @defcolumn,
    @deftable,
    @defgroup,
    insert_column!,
    insert_table!,
    group_table!,
    column_creator,
    table_creator,
    table_grouper,
    covid_19_database,

    PARAMS,
    NPARAMS,
    VARIABLES,
    NVARIABLES,
    ParamsNamedTuple,
    VariablesNamedTuple,
    model_ode,
    model_ode_function,
    pack_variables,
    unpack_variables,
    pack_parameters,
    unpack_parameters,
    model_parameters,
    model_initial_values,
    model_problem,
    model_solution,
    model_decode_solution,
    model_plot,
    model_loss_function

const DATA_SOURCES = Dict{Symbol, Type}()

include("utils.jl")
include("datasources.jl")
include("database.jl")
include("model.jl")

end # module
