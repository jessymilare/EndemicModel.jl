# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

module EndemicModel

import Base.Iterators
import CSV
import JSON
import Measures: mm
import Plots
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
using Gettext
using GZip
using Interpolations
using Optim
using StatsBase
using Unitful

export DataDict,
    DataFrameDict,
    AbstractDataDict,
    AbstractDataSource,
    DataSourceDesignator,
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
    SEIRVariables,
    SEIRDerivatives,
    SEIRParameters,
    SEIRModel,
    SEIRModel!,
    pack_vars,
    pack_params,
    unpack_vars,
    unpack_params,
    SEIR_ODE_fun,
    SEIR_ODEProblem,
    plot,
    model_problem,
    model_solution,
    to_dataframe,
    to_dataframe!,
    model_validate,
    model_step,
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
    realdata,
    realdata!,
    phuber_loss,
    diff_phuber_loss,
    l2_loss,
    diff_l2_loss,
    model_loss,
    optimize_parameters,
    optimize_parameters!,
    optimize_variables,
    optimize_variables!,
    estimate_μ,
    estimate_α,
    estimate_γ,
    estimate_β,
    estimate_exposed!,
    covid19_database,
    covid19

const DATA_SOURCES = Dict{Symbol,Type}()

include("utils.jl")
include("options.jl")
include("source.jl")
include("sources/general.jl")
include("sources/covid19.jl")
include("database.jl")
include("subdatabase.jl")
include("database_files.jl")
include("covid19.jl")
include("model.jl")
include("plot.jl")
include("estimate.jl")

end # module
