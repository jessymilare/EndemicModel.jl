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
function modeldataof(model::AbstractEndemicModel)
    model.modeldata
end
function modeldataof!(model::AbstractEndemicModel, value::AbstractDataFrame)
    model.modeldata = value
end

dataof(model::AbstractEndemicModel) = model.data
dataof!(model::AbstractEndemicModel, value) = model.data = value

const SEIR_VARS = (:S, :E, :I, :R)
const SEIR_DERIV = (:dS, :dE, :dI, :dR)
const SEIR_PARAMS = (:M, :β, :γ, :α, :μ)
const SEIRVars = NamedTuple{SEIR_VARS, NTuple{4, Vector{Float64}}}
const SEIRDeriv = NamedTuple{SEIR_DERIV, NTuple{4, Vector{Float64}}}
const SEIRParams = NamedTuple{SEIR_PARAMS, NTuple{5, Vector{Float64}}}

SEIRVars(S, E, I, R) = SEIRVars((S, E, I, R))
SEIRParams(M, β, γ, α, μ) = SEIRParams((M, β, γ, α, μ))

function SEIRDeriv(vars::SEIRVars, params::SEIRParams)
    (S, E, I, R) = vars
    (M, β, γ, α, μ) = params

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
dS = - β * (I ./ sum(M)) .* S
dE = β * (I ./ sum(M)) .* S - γ .* E
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
    StoE = β .* sum(I ./ pop) .* S
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
    vars::SEIRVars
    deriv::SEIRDeriv
    params::SEIRParams
    ngroups::Int
    groupnames::Vector{Symbol}
    data::Union{Nothing, DataFrame}
    modeldata::Union{Nothing, DataFrame}
    info::Dict{Symbol, Any}

    function SEIRModel(
        vars::SEIRVars,
        params::SEIRParams,
        groupnames::Union{Nothing, Vector{Symbol}} = nothing;
        kwargs...,
    )
        ngroups = length(vars[1])
        if groupnames == nothing
            if ngroups == 1
                groupnames = [:all]
            else
                @error "Parameter `groupnames` must be provided."
            end
        end

        deriv = SEIRDeriv(vars, params)
        info = Dict(kwargs...)

        new(vars, deriv, params, ngroups, groupnames, nothing, nothing, info)
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
    μ::AbstractVector{Float64},
    groupnames::Union{Nothing, Vector{Symbol}} = nothing;
    kwargs...,
)
    SEIRModel(SEIRVars(S, E, I, R), SEIRParams(M, β, γ, α, μ), groupnames; kwargs...)
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

function pack_params(M, β, γ, α, μ)
    args = filter(!isnothing, (M, β, γ, α, μ))
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
        μ = :μ ∈ packed_params ? tup[:μ] : nothing
        pack_params(M, β, γ, α, μ)
    end
end

function unpack_vars(X::Vector{Float64})
    nv = 4
    Xlen = length(X)
    S = view(X, 1:nv:Xlen)
    E = view(X, 2:nv:Xlen)
    I = view(X, 3:nv:Xlen)
    R = view(X, 4:nv:Xlen)
    SEIRVars(S, E, I, R)
end

function unpack_params(
    P::Vector{Float64},
    ngroups::Int;
    packed_params = SEIR_PARAMS,
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
    SEIRParams(M, β, γ, α, μ)
end

function SEIR_ODE_fun(X, P, t)
    (S, E, I, R) = unpack_vars(X)
    n = length(S)
    (M, β, γ, α, μ) = unpack_params(P, n)

    (dS, dE, dI, dR) = SEIRDeriv(SEIRVars(S, E, I, R), SEIRParams(M, β, γ, α, μ))
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

function to_dataframe(problem::ODEProblem, model::SEIRModel; kwargs...)
    solution = solve(problem, Euler(); dt = 1.0)
    (M, β, γ, α, μ) = paramsof(model)
    days = solution.t
    initial_date = get(model.info, :initial_date, model.data.date[end])
    date = initial_date .+ Day.(days)
    n = length(days)
    ng = model.ngroups
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
        recovered[i] = sum(R .* (1.0 .- μ))
        deaths[i] = sum(R .* μ)
        confirmed[i] = active[i] + recovered[i] + deaths[i]
        diff_exposed[i] = sum(dE)
        diff_active[i] = sum(dI)
        diff_recovered[i] = sum(dR .* (1.0 .- μ))
        diff_deaths[i] = sum(dR .* μ)
    end
    DataFrame(
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
end

function to_dataframe!(model::SEIRModel; kwargs...)
    modeldataof!(model, to_dataframe(model; kwargs...))
end

function to_dataframe(model::SEIRModel; kwargs...)
    to_dataframe(model_problem(model; kwargs...), model; kwargs...)
end

function model_plot(
    df::DataFrame;
    plot_columns = option(:plot_columns),
    title = summary(df),
    yfactor = 1e-6,
    ylabel = yfactor == 1e-6 ? "People (millions)" : "People",
    minimum_plot_factor = option(:minimum_plot_factor),
    kwargs...,
)
    colnames = intersect(plot_columns, Symbol.(names(df)))
    numpeople = if hasproperty(df, :estimated_population)
        df.estimated_population[1]
    else
        df.confirmed[end]
    end
    n = findlast(df.active .>= numpeople * minimum_plot_factor)
    df = df[1:n, :]
    xlabel = :date
    y1 = colnames[1]
    win = plot(
        df[!, xlabel],
        df[!, y1] .* yfactor,
        xlabel = string(xlabel),
        label = string(y1),
        ylabel = ylabel,
        title = title,
        legend = :right,
    )
    for yn ∈ colnames[2:end]
        ynstr = string(yn)
        plot!(
            win,
            df[!, xlabel],
            df[!, yn] .* yfactor;
            xlabel = string(xlabel),
            label = ynstr,
            ylabel = ylabel,
        )
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
