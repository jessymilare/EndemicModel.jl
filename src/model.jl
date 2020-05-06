abstract type AbstractEndemicModel end

function varsof(model::AbstractEndemicModel)
    model.vars
end
function varsof!(model::AbstractEndemicModel, value)
    model.vars = value
end
function derivof(model::AbstractEndemicModel)
    model.deriv
end
function derivof!(model::AbstractEndemicModel, value)
    model.deriv = value
end
function paramsof(model::AbstractEndemicModel)
    model.params
end
function paramsof!(model::AbstractEndemicModel, value)
    model.params = value
end
function modeldata_of(model::AbstractEndemicModel)
    model.modeldata
end
function modeldata_of!(model::AbstractEndemicModel, value::AbstractDataFrame)
    model.modeldata = value
end

dataof(model::AbstractEndemicModel) = model.data

const SEIR_VARS = (:S, :E, :I, :R)
const SEIR_DERIV = (:dS, :dE, :dI, :dR)
const SEIR_PARAMS = (:M, :β, :γ, :α)
const SEIRVars = NamedTuple{SEIR_VARS, NTuple{4, Vector{Float64}}}
const SEIRDeriv = NamedTuple{SEIR_DERIV, NTuple{4, Vector{Float64}}}
const SEIRParams = NamedTuple{SEIR_PARAMS, NTuple{4, Vector{Float64}}}

SEIRVars(S, E, I, R) = SEIRVars((S, E, I, R))
SEIRParams(M, β, γ, α) = SEIRParams((M, β, γ, α))

function SEIRDeriv(vars::SEIRVars, params::SEIRParams)
    (S, E, I, R) = vars
    (M, β, γ, α) = params

    SEIRDeriv(S, E, I, R, M, β, γ, α)
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
dS = - β * ((E .+ I) ./ sum(M)) .* S
dE = β * ((E .+ I) ./ sum(M)) .* S - γ .* E
dI = γ .* E - α .* I
dR = α .* I
"""
function SEIRDeriv(
    S::AbstractVector{Float64},
    E::AbstractVector{Float64},
    I::AbstractVector{Float64},
    R::AbstractVector{Float64},
    M::AbstractVector{Float64},
    β::AbstractVector{Float64},
    γ::AbstractVector{Float64},
    α::AbstractVector{Float64},
)
    pop = sum(M)
    StoE = β .* sum((E .+ I) ./ pop) .* S
    EtoI = γ .* E
    ItoR = α .* I

    dS = max.(-StoE, -S)
    dE = max.(StoE - EtoI, -E)
    dI = max.(EtoI - ItoR, -I)
    dR = max.(ItoR, -R)

    (; dS = dS, dE = dE, dI = dI, dR = dR)
end

initial_M(kind::Symbol; kwargs...) = initial_M(Val(kind); kwargs...)

function initial_M(kind::Val{:nogroup}; estimated_population = 1e6, kwargs...)
    Float64[estimated_population]
end

function initial_M(
    kind::Val{:lethal_group};
    estimated_population = 1e6,
    μ = μ_covid_19,
    kwargs...,
)
    Float64(estimated_population) .* Float64[1.0 - μ, μ]
end

function initial_M(
    kind::Val{:serious_lethal_groups};
    estimated_population = 1e6,
    μ = μ_covid_19,
    serious_rate = 0.1,
    kwargs...,
)
    Float64(estimated_population) .* [1.0 - μ - serious_rate, serious_rate, μ]
end

function initial_M(
    kind::Val{:lethal_conf_groups};
    estimated_population = 1e6,
    μ = μ_covid_19,
    kwargs...,
)
    Float64(estimated_population) .* Float64[1.0 - μ, μ, 0.0, 0.0]
end

"""
Represent a SEIR differential equation.

S, E, I, R must be vectors of the same length `n`.
M, γ must be vectors of length `n`.
β must be an array of length `n x n`.

Variables:
S = susceptible
E = exposed (infected, but not tested and confirmed)
I = active (confirmed)
R = recovered / deaths

Parameters:
M = total susceptible population
β = transmission rate
γ = evolution rate, i.e., E -> I rate
α = recovery/death rate, i.e., Ic -> R rate

Equations:
M = S + E + I + R
dM = 0
dS = - β * (I ./ M) .* S
dE = β * (I ./ M) .* S - γ .* E
dI = γ .* E - α .* I
dR = α .* I
"""
mutable struct SEIRModel <: AbstractEndemicModel
    vars::SEIRVars
    deriv::SEIRDeriv
    params::SEIRParams
    ngroups::Int
    groupnames::Vector{Symbol}
    lethal_groups::Vector{Int}
    unconfirmed_groups::Vector{Int}
    data::Union{Nothing, DataFrame}
    modeldata::Union{Nothing, DataFrame}

    function SEIRModel(
        vars::SEIRVars,
        params::SEIRParams,
        groupnames::Union{Nothing, Vector{Symbol}} = nothing,
    )
        ngroups = length(vars[1])
        if groupnames == nothing
            if ngroups == 1
                groupnames = [:all]
                lethal_groups = Int[]
                unconfirmed_groups = Int[]
            elseif ngroups == 2
                groupnames = [:non_lethal, :lethal]
                lethal_groups = Int[2]
                unconfirmed_groups = Int[]
            elseif ngroups == 3
                groupnames = [:mild, :serious, :lethal]
                lethal_groups = Int[3]
                unconfirmed_groups = Int[]
            else
                @error "Parameter `groupnames` must be provided."
            end
        else
            lethal_groups = Int[]
            unconfirmed_groups = Int[]
            for i ∈ 1:ngroups
                if string(groupnames[i])[1:6] == "lethal"
                    push!(lethal_groups, i)
                end
            end
        end

        deriv = SEIRDeriv(vars, params)

        new(
            vars,
            deriv,
            params,
            ngroups,
            groupnames,
            lethal_groups,
            unconfirmed_groups,
            nothing,
            nothing,
        )
    end
end

function SEIRModel(
    S::AbstractVector{Float64},
    E::AbstractVector{Float64},
    I::AbstractVector{Float64},
    R::AbstractVector{Float64},
    M::AbstractVector{Float64},
    β::AbstractVector{Float64},
    γ::AbstractVector{Float64},
    α::AbstractVector{Float64},
    groupnames::Union{Nothing, Vector{Symbol}} = nothing,
)
    SEIRModel(
        (; S = S, E = E, I = I, R = R),
        (; M = M, β = β, γ = γ, α = α),
        groupnames,
    )
end

function step(model::SEIRModel)
    (S, E, I, R) = varsof(model)
    (dS, dE, dI, dR) = derivof(model)
    vars = SEIRVars(S + dS, E + dE, I + dI, R + dR)
    params = paramsof(model)
    SEIRModel(vars, params)
end

function step!(model::SEIRModel)
    (S, E, I, R) = varsof(model)
    (dS, dE, dI, dR) = derivof(model)
    vars = SEIRVars(S + dS, E + dE, I + dI, R + dR)
    params = paramsof(model)
    deriv = SEIRDeriv(vars, params)
    varsof!(model, vars)
    derivof!(model, deriv)
    model
end

function pack_vars(S, E, I, R)
    vectors = zip(S, E, I, R)
    result = Float64[]
    for vec ∈ vectors
        append!(result, vec)
    end
    result
end

function pack_vars(tup::Union{Tuple, NamedTuple})
    pack_vars(tup...)
end

function pack_params(M, β, γ, α)
    args = filter(!isnothing, (M, β, γ, α))
    vectors = zip(args...)
    result = Float64[]
    for v ∈ vectors
        append!(result, v)
    end
    result
end

function pack_params(tup::Tuple)
    pack_params(tup...)
end

function pack_params(tup::NamedTuple; packed_params = keys(tup))
    if packed_params == SEIR_PARAMS
        pack_params(tup...)
    else
        M = :M ∈ packed_params ? tup[:M] : nothing
        β = :β ∈ packed_params ? tup[:β] : nothing
        γ = :γ ∈ packed_params ? tup[:γ] : nothing
        α = :α ∈ packed_params ? tup[:α] : nothing
        pack_params(M, β, γ, α)
    end
end

function unpack_vars(X::Vector{Float64})
    nv = 4
    Xlen = length(X)
    S = view(X, 1:nv:Xlen)
    E = view(X, 2:nv:Xlen)
    I = view(X, 3:nv:Xlen)
    R = view(X, 4:nv:Xlen)
    (; S = S, E = E, I = I, R = R)
end

function unpack_params(
    P::Vector{Float64},
    ngroups::Int;
    packed_params = SEIR_PARAMS,
    default_params = (),
)
    vecvars = [:M, :β, :γ, :α] ∩ packed_params
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
    (; M = M, β = β, γ = γ, α = α)
end

function SEIR_ODE_fun(X, P, t)
    (S, E, I, R) = unpack_vars(X)
    n = length(S)
    (M, β, γ, α) = unpack_params(P, n)

    (dS, dE, dI, dR) = SEIRDeriv(SEIRVars(S, E, I, R), SEIRParams(M, β, γ, α))
    pack_vars(dS, dE, dI, dR)
end

function SEIR_ODEProblem(
    inivars::SEIRVars,
    params::SEIRParams;
    maxtime = 180.0,
    kwargs...,
)
    ODEProblem(SEIR_ODE_fun, pack_vars(inivars), (0.0, maxtime), pack_params(params))
end

function model_problem(model::SEIRModel; kwargs...)
    SEIR_ODEProblem(varsof(model), paramsof(model); kwargs...)
end

function model_solution(model::SEIRModel; kwargs...)
    solve(model_problem(model; kwargs...), Euler(); dt = 1.0)
end

function to_dataframe(problem::ODEProblem, model::SEIRModel)
    solution = solve(problem, Euler(); dt = 1.0)
    (M, β, γ, α) = paramsof(model)
    days = solution.t
    date = model.data.date[end] .+ Day.(days)
    n = length(days)
    ng = model.ngroups
    lethal = model.lethal_groups
    exposed = zeros(Float64, n)
    active = zeros(Float64, n)
    recovered = zeros(Float64, n)
    deaths = zeros(Float64, n)
    confirmed = zeros(Float64, n)
    diff_exposed = zeros(Float64, n)
    diff_active = zeros(Float64, n)
    diff_recovered = zeros(Float64, n)
    diff_deaths = zeros(Float64, n)
    for i ∈ 1:n
        (S, E, I, R) = unpack_vars(solution[i])
        (dS, dE, dI, dR) = SEIRDeriv(S, E, I, R, M, β, γ, α)
        exposed[i] = sum(E)
        active[i] = sum(I)
        recovered[i] = sum(R[setdiff(1:ng, lethal)])
        deaths[i] = sum(R[lethal])
        confirmed[i] = active[i] + recovered[i] + deaths[i]
        diff_exposed[i] = sum(dE)
        diff_active[i] = sum(dI)
        diff_recovered[i] = sum(dR[setdiff(1:ng, lethal)])
        diff_deaths[i] = sum(dR[lethal])
    end
    df = DataFrame(
        date = date,
        confirmed = confirmed,
        exposed = exposed,
        active = active,
        recovered = recovered,
        deaths = deaths,
        diff_exposed = diff_exposed,
        diff_active = diff_active,
        diff_recovered = diff_recovered,
        diff_deaths = diff_deaths,
    )

    df
end

function to_dataframe!(model::SEIRModel; kwargs...)
    model.modeldata = to_dataframe(model_problem(model; kwargs...), model)
end

function to_dataframe(model::SEIRModel; kwargs...)
    to_dataframe(model_problem(model; kwargs...), model)
end

const PLOT_COLUMNS = [:date, :confirmed, :exposed, :active, :recovered, :deaths]

function model_plot(
    df::DataFrame,
    colnames = intersect(names(df), PLOT_COLUMNS);
    title = summary(df),
    minimum_plot_factor = option(:minimum_plot_factor),
    kwargs...,
)
    numpeople = if hasproperty(df, :estimated_population)
        df.estimated_population[1]
    else
        df.confirmed[end]
    end
    n = findlast(df.active .>= numpeople * minimum_plot_factor)
    df = df[1:n, :]
    xlabel = colnames[1]
    y1 = colnames[2]
    win = plot(
        df[!, xlabel],
        df[!, y1],
        xlabel = string(xlabel),
        label = string(y1),
        title = title,
        legend = :topleft,
    )
    for yn ∈ colnames[3:end]
        ynstr = string(yn)
        plot!(win, df[!, xlabel], df[!, yn]; xlabel = string(xlabel), label = ynstr)
    end
    gui(win)
end

function model_plot(model::SEIRModel; kwargs...)
    df = to_dataframe!(model; kwargs...)
    model_plot(df; kwargs...)
end

function model_plot(problem::ODEProblem; kwargs...)
    df = to_dataframe(problem; kwargs...)
    model_plot(df; kwargs...)
end
