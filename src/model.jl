# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

abstract type AbstractEndemicModel end

variables(model::AbstractEndemicModel) = model.variables
variables!(model::AbstractEndemicModel, value) = model.variables = value

derivatives(model::AbstractEndemicModel) = model.derivatives
derivatives!(model::AbstractEndemicModel, value) = model.derivatives = value

parameters(model::AbstractEndemicModel) = model.parameters
parameters!(model::AbstractEndemicModel, value) = model.parameters = value

modeldata(model::AbstractEndemicModel) = model.modeldata
modeldata!(model::AbstractEndemicModel, value::AbstractDataFrame) = model.modeldata = value

modeldata!(model::AbstractEndemicModel; kwargs...) =
    modeldata!(model::AbstractEndemicModel, to_dataframe(model; kwargs...))

default_kwargs(model::AbstractEndemicModel) = model.kwargs
default_kwargs!(model::AbstractEndemicModel, value) = model.kwargs = value

default_kwarg(model::AbstractEndemicModel, key::Symbol) = model.kwargs[key]
default_kwarg(model::AbstractEndemicModel, key::Symbol, default) =
    get(model.kwargs, key, default)
default_kwarg!(model::AbstractEndemicModel, key::Symbol, value) = model.kwargs[key] = value

realdata(model::AbstractEndemicModel) = model.realdata
realdata!(model::AbstractEndemicModel, value) = model.realdata = value

groupnames(model::AbstractEndemicModel) = model.groupnames
groupnames!(model::AbstractEndemicModel, value) = model.groupnames = value

ngroups(model::AbstractEndemicModel) = model.ngroups
ngroups!(model::AbstractEndemicModel, value) = model.ngroups = value

export_data(model::AbstractEndemicModel; kwargs...) =
    export_data(modeldata(model); kwargs...)

export_data(source, model::AbstractEndemicModel; kwargs...) =
    export_data(source, modeldata(model); kwargs...)

const SEIR_VARIABLES = (:S, :E, :I, :R)
const SEIR_DERIVATIVES = (:dS, :dE, :dI, :dR)
const SEIR_PARAMETERS = (:M, :β, :γ, :α, :μ)
const SEIRVariables = NamedTuple{SEIR_VARIABLES,NTuple{4,Vector{Float}}}
const SEIRDerivatives = NamedTuple{SEIR_DERIVATIVES,NTuple{4,Vector{Float}}}
const SEIRParameters = NamedTuple{SEIR_PARAMETERS,NTuple{5,Vector{Float}}}

SEIRVariables(S, E, I, R) = SEIRVariables((S, E, I, R))
SEIRDerivatives(dS, dE, dI, dR) = SEIRDerivatives((dS, dE, dI, dR))
SEIRParameters(M, β, γ, α, μ) = SEIRParameters((M, β, γ, α, μ))

function _SEIR_derivative(vars::SEIRVariables, params::SEIRParameters)
    (S, E, I, R) = vars
    (M, β, γ, α, μ) = params

    _SEIR_derivative(S, E, I, R, M, β, γ, α)
end

"""
Compute the derivative of a SEIR-type differential equation.
In this system, `exposed` people are also considered infections, yet not tested confirmed.
`S`, `E`, `I`, `R` must be vectors of the same length `n`.
`M`, `γ` must be vectors of length `n`.
`β` must be Float64.

Equations:
M = S + E + I + R
dM = 0
dS = - β * (I ./ sum(M)) .* S
dE = β * (I ./ sum(M)) .* S - γ .* E
dI = γ .* E - α .* I
dR = α .* I
"""
function _SEIR_derivative(
    S::AbstractVector{Float},
    E::AbstractVector{Float},
    I::AbstractVector{Float},
    R::AbstractVector{Float},
    M::AbstractVector{Float},
    β::AbstractVector{Float},
    γ::AbstractVector{Float},
    α::AbstractVector{Float};
    backward_limit::Bool = false,
)
    pop = sum(M)
    StoE = β .* sum(I ./ pop) .* S
    EtoI = γ .* E
    ItoR = α .* I
    if backward_limit
        StoE = min.(E, StoE)
        EtoI = min.(I, EtoI)
        ItoR = min.(R, ItoR)
    end

    dS = -StoE
    dE = StoE - EtoI
    dI = EtoI - ItoR
    dR = ItoR

    # @debug "Computed SEIR derivatives" (S, E, I, R) (dS, dE, dI, dR)

    (; dS = dS, dE = dE, dI = dI, dR = dR)
end

initial_M(kind::Symbol; kwargs...) = initial_M(Val(kind); kwargs...)

function initial_M(kind::Val{:nogroup}; estimated_population = 1e6, kwargs...)
    Float[estimated_population]
end

"""
Represent a SEIR differential equation.

S, E, I, R must be vectors of the same length `n`.
M, β, γ, α, μ must be vectors of length `n`.

Variables:
S = susceptible
E = exposed (infected, but not infeccious)
I = active (infeccious, confirmed)
R = recovered / deaths

Parameters:
M = total susceptible population
β = transmission rate
γ = evolution rate, i.e., E -> I rate
α = recovery/death rate, i.e., Ic -> R rate
μ = death rate

Equations:
M = S + E + I + R
dM = 0
dS = - β * (I ./ M) .* S
dE = β * (I ./ M) .* S - γ .* E
dI = γ .* E - α .* I
dR = α .* I
"""
mutable struct SEIRModel <: AbstractEndemicModel
    variables::SEIRVariables
    parameters::SEIRParameters
    derivatives::SEIRDerivatives
    ngroups::Int
    groupnames::Vector{Symbol}
    realdata::Union{Nothing,AbstractDataFrame}
    modeldata::Union{Nothing,AbstractDataFrame}
    kwargs::Dict{Symbol,Any}

    function SEIRModel(
        vars::SEIRVariables,
        params::SEIRParameters;
        ngroups::Int = length(vars[1]),
        groupnames = nothing,
        realdata::Union{Nothing,AbstractDataFrame} = nothing,
        modeldata::Union{Nothing,AbstractDataFrame} = nothing,
        kwargs...,
    )
        groupnames =
            something(groupnames, ngroups == 1 ? [:all] : "group_" .* string.(1:ngroups))
        groupnames = collect(Symbol.(groupnames))
        deriv = _SEIR_derivative(vars, params)
        info = Dict{Symbol,Any}(kwargs)

        model = new(vars, params, deriv, ngroups, groupnames, realdata, modeldata, info)
        if isnothing(modeldata)
            modeldata!(model, to_dataframe(model; kwargs...))
        end
        model
    end
end

function Base.copy(model::SEIRModel)
    SEIRModel(
        variables(model),
        parameters(model);
        ngroups = ngroups(model),
        groupnames = groupnames(model),
        realdata = realdata(model),
        modeldata = modeldata(model),
        default_kwargs(model)...,
    )
end

SEIRModel(model::SEIRModel) = model

function SEIRModel(
    S::AbstractVector{Float},
    E::AbstractVector{Float},
    I::AbstractVector{Float},
    R::AbstractVector{Float},
    M::AbstractVector{Float},
    β::AbstractVector{Float},
    γ::AbstractVector{Float},
    α::AbstractVector{Float},
    μ::AbstractVector{Float};
    kwargs...,
)
    SEIRModel(SEIRVariables(S, E, I, R), SEIRParameters(M, β, γ, α, μ); kwargs...)
end

function model_step(model::SEIRModel, n::Int = 1; initial_date = nothing, kwargs...)
    (S, E, I, R) = variables(model)
    (dS, dE, dI, dR) = derivatives(model)
    params = parameters(model)
    (M, β, γ, α, μ) = params
    initial_date = something(
        initial_date,
        default_kwarg(model, :initial_date, today() - Day(15)) + Day(n),
    )
    data = realdata(model)
    if n > 0
        while n > 0
            for i ∈ 1:4
                (S, E, I, R) = (S + 0.25dS, E + 0.25dE, I + 0.25dI, R + 0.25dR)
                (dS, dE, dI, dR) = _SEIR_derivative(S, E, I, R, M, β, γ, α)
            end
            n -= 1
        end
    else
        z = zeros(Float, length(S))
        if !isnothing(data)
            idx = Dates.days(initial_date - data.date[1]) + 1
            if idx >= 1
                n = 0 # skip loop below
            else
                n, idx = idx - 1, 1
            end
            E = [data.exposed[idx]]
            I = [data.active[idx]]
            R = [data.closed[idx]]
            S = [data.estimated_population[idx] - E[1] - I[1] - R[1]]
        end
        while n < 0
            for i ∈ 1:4
                (S, E, I, R) = (S - 0.25dS, E - 0.25dE, I - 0.25dI, R - 0.25dR)
                (dS, dE, dI, dR) =
                    _SEIR_derivative(S, E, I, R, M, β, γ, α; backward_limit = true)
            end
            n += 1
        end

        # optimize_variables!(model; kwargs...)
        modeldata!(model, to_dataframe(model; kwargs...))
    end
    SEIRModel(
        SEIRVariables(S, E, I, R),
        params;
        initial_date = initial_date,
        realdata = data,
        kwargs...,
    )
end

function pack_vars(S, E, I, R)
    vectors = zip(S, E, I, R)
    result = Float[]
    for vec ∈ vectors
        append!(result, vec)
    end
    result
end

function pack_vars(tup::Union{Tuple,NamedTuple})
    pack_vars(tup...)
end

function pack_params(M, β, γ, α, μ)
    args = filter(!isnothing, (M, β, γ, α, μ))
    vectors = zip(args...)
    result = Float[]
    for v ∈ vectors
        append!(result, v)
    end
    result
end

function pack_params(tup::Tuple)
    pack_params(tup...)
end

function pack_params(tup::NamedTuple; packed_params = keys(tup))
    if packed_params == SEIR_PARAMETERS
        pack_params(tup...)
    else
        M = :M ∈ packed_params ? tup[:M] : nothing
        β = :β ∈ packed_params ? tup[:β] : nothing
        γ = :γ ∈ packed_params ? tup[:γ] : nothing
        α = :α ∈ packed_params ? tup[:α] : nothing
        μ = :μ ∈ packed_params ? tup[:μ] : nothing
        pack_params(M, β, γ, α, μ)
    end
end

function unpack_vars(X::Vector{Float})
    nv = 4
    Xlen = length(X)
    S = view(X, 1:nv:Xlen)
    E = view(X, 2:nv:Xlen)
    I = view(X, 3:nv:Xlen)
    R = view(X, 4:nv:Xlen)
    SEIRVariables(S, E, I, R)
end

function unpack_params(
    P::Vector{Float},
    ngroups::Int;
    packed_params = SEIR_PARAMETERS,
    default_params = (),
)
    vecvars = [:M, :β, :γ, :α, :μ] ∩ packed_params
    Plen = length(vecvars) * ngroups
    index = 1
    np = length(vecvars)
    if :M ∈ vecvars
        M = view(P, index:np:Plen)
        index += 1
    else
        M = default_params[:M]
    end
    if :β ∈ vecvars
        β = view(P, index:np:Plen)
        index += 1
    else
        β = default_params[:M]
    end
    if :γ ∈ vecvars
        γ = view(P, index:np:Plen)
        index += 1
    else
        γ = default_params[:γ]
    end
    if :α ∈ vecvars
        α = view(P, index:np:Plen)
        index += 1
    else
        α = default_params[:α]
    end
    if :μ ∈ vecvars
        μ = view(P, index:np:Plen)
        index += 1
    else
        μ = default_params[:μ]
    end
    SEIRParameters(M, β, γ, α, μ)
end

function SEIR_ODE_fun(X, P, t)
    (S, E, I, R) = unpack_vars(X)
    n = length(S)
    (M, β, γ, α, μ) = unpack_params(P, n)

    (dS, dE, dI, dR) =
        _SEIR_derivative(SEIRVariables(S, E, I, R), SEIRParameters(M, β, γ, α, μ))
    pack_vars(dS, dE, dI, dR)
end

function SEIR_ODEProblem(
    inivars::SEIRVariables,
    params::SEIRParameters;
    maxtime = 360,
    kwargs...,
)
    ODEProblem(SEIR_ODE_fun, pack_vars(inivars), (0.0, Float(maxtime)), pack_params(params))
end

function model_problem(model::SEIRModel; kwargs...)
    SEIR_ODEProblem(variables(model), parameters(model); kwargs...)
end

function model_solution(model::SEIRModel; dt = 1.0, kwargs...)
    solve(model_problem(model; kwargs...), Euler(); dt = dt)
end

function to_dataframe(model::SEIRModel; diff_columns::Bool = true, kwargs...)
    solution = model_solution(model; dt = 0.25, kwargs...)
    (M, β, γ, α, μ) = parameters(model)
    days = solution.t[1:4:end]
    solution = solution[1:4:end]
    initial_date = default_kwarg(model, :initial_date, today() - Day(15))
    date = initial_date .+ Day.(days)
    n = length(days)
    ng = model.ngroups
    exposed = zeros(Int, n)
    active = zeros(Int, n)
    recovered = zeros(Int, n)
    deaths = zeros(Int, n)
    infected = zeros(Int, n)
    if diff_columns
        diff_exposed = zeros(Int, n)
        diff_active = zeros(Int, n)
        diff_recovered = zeros(Int, n)
        diff_deaths = zeros(Int, n)
    end

    _round(x) = isnan(x) ? 0 : round(Int, max(0, min(1e12, x)))

    for i ∈ 1:n
        (S, E, I, R) = unpack_vars(solution[i])
        (dS, dE, dI, dR) = _SEIR_derivative(S, E, I, R, M, β, γ, α)
        exposed[i] = sum(E) |> _round
        active[i] = sum(I) |> _round
        recovered[i] = sum(R .* (1.0 .- μ)) |> _round
        deaths[i] = sum(R .* μ) |> _round
        infected[i] = active[i] + recovered[i] + deaths[i] |> _round
        if diff_columns
            diff_exposed[i] = sum(dE) |> _round
            diff_active[i] = sum(dI) |> _round
            diff_recovered[i] = sum(dR .* (1.0 .- μ)) |> _round
            diff_deaths[i] = sum(dR .* μ) |> _round
        end
    end
    dcols = if diff_columns
        [
            :diff_exposed => diff_exposed,
            :diff_active => diff_active,
            :diff_recovered => diff_recovered,
            :diff_deaths => diff_deaths,
        ]
    else
        []
    end
    DataFrame(
        :date => date,
        :infected => infected,
        :exposed => exposed,
        :active => active,
        :recovered => recovered,
        :deaths => deaths,
        dcols...;
        copycols = false,
    )
end
