module EndemicModel

import Base.Iterators
import CSV
import JSON
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
    to_json,

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

    SEIRModel,
    pack_vars,
    pack_params,
    unpack_vars,
    unpack_params,
    SEIR_ODE_fun,
    SEIR_ODEProblem,
    model_plot,
    model_problem,
    model_solution,
    to_dataframe,
    to_dataframe!,

    dataof,
    dataof!,
    paramsof,
    paramsof!,

    phuber_loss,
    diff_phuber_loss,
    l2_loss,
    diff_l2_loss,
    model_loss,
    optimize_params,
    estimate_μ,
    estimate_α,
    estimate_γ,
    estimate_β

const DATA_SOURCES = Dict{Symbol, Type}()

include("utils.jl")
include("datasources.jl")
include("database.jl")
include("covid19.jl")
include("model.jl")
include("estimate.jl")

end # module
