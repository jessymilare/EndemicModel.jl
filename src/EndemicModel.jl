# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

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
using StatsBase
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

    SEIR_PARAMS,
    SEIR_DERIV,
    SEIR_VARS,
    SEIRVars,
    SEIRDeriv,
    SEIRParams,
    SEIRModel,
    SEIRModel!,
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

    datadict,
    datadict!,
    modeldict,
    modeldict!,
    paramdict,
    paramdict!,

    variables,
    variables!,
    derivatives,
    derivatives!,
    parameters,
    parameters!,
    modeldata,
    modeldata!,

    phuber_loss,
    diff_phuber_loss,
    l2_loss,
    diff_l2_loss,
    model_loss,
    optimize_params,
    estimate_μ,
    estimate_α,
    estimate_γ,
    estimate_β,
    estimate_exposed!,

    covid_19_database,
    covid_19

const DATA_SOURCES = Dict{Symbol, Type}()

include("utils.jl")
include("options.jl")
include("datasources.jl")
include("database.jl")
include("database_files.jl")
include("covid19.jl")
include("model.jl")
include("estimate.jl")

end # module
