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

dataof(model::AbstractEndemicModel) = model.data

const SEIR_VARS = (:S, :E, :I, :R)
const SEIR_DERIV = (:dS, :dE, :dI, :dR)
const SEIR_PARAMS = (:M, :β, :γ, :α)
const SEIRVars = NamedTuple{SEIR_VARS, NTuple{4, Vector{Float64}}}
const SEIRDeriv = NamedTuple{SEIR_DERIV, NTuple{4, Vector{Float64}}}
const SEIRParams = NamedTuple{SEIR_PARAMS, NTuple{4, Vector{Float64}}}

SEIRVars(S, E, I, R) = (SEIRVars((S, E, I, R)))
SEIRParams(M, β, γ, α) = (SEIRParams((M, β, γ, α)))
SEIRDeriv(dS, dE, dI, dR) = SEIRDeriv((dS, dE, dI, dR))

function SEIRDeriv(vars::SEIRVars, params::SEIRParams)
    (S, E, I, R) = vars
    (M, β, γ, α) = params
    StoE = sum(β .* (I ./ M)) * S
    EtoI = γ .* E
    ItoR = α .* I

    dS = -StoE
    dE = StoE - EtoI
    dI = EtoI - ItoR
    dR = ItoR

    (; dS = dS, dE = dE, dI = dI, dR = dR)
end

"""
Compute the derivative of a SEIR-type differential equation.
S, E, I, R must be vectors of the same length `n`.
M, γ must be vectors of length `n`.
β must be an array of length `n x n`.

Equations:
M = S + E + I + R
dM = 0
dS = - β * (I ./ M) .* S
dE = β * (I ./ M) .* S - γ .* E
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
    SEIRDeriv((; S = S, E = E, I = I, R = R), (; M = M, β = β, γ = γ, α = α))
end

"""
Represent a SEIR differential equation.

S, E, I, R must be vectors of the same length `n`.
M, γ must be vectors of length `n`.
β must be an array of length `n x n`.

Variables:
S = susceptible
E = exposed (non-active during incubation time)
I =  active
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
            elseif ngroups == 4
                groupnames = [
                    :non_lethal_unconfirmed,
                    :lethal_unconfirmed,
                    :non_lethal_confirmed,
                    :lethal_confirmed,
                ]
                lethal_groups = Int[4]
                unconfirmed_groups = Int[1, 2]
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
                if string(groupnames[i])[(end - 11):end] == "unconfirmed"
                    push!(unconfirmed_groups, i)
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

function SEIRModel(data::DataFrame)
    numpeople = data.estimated_population[1]
    deaths = data.deaths[end]
    recovered = data.recovered[end]
    μ = deaths / (recovered + deaths)
    M = numpeople * [1 - μ, μ]
    I = data.active[1] * [1 - μ, μ]
    S = M - I
    E = Float64[0.0, 0.0]
    R = Float64[data.recovered[1], data.deaths[1]]
    β = [0.15, 0.15]
    γ = [1.0 / 3.0, 1.0 / 3.0]
    α = [1.0 / 17.0, 1.0 / 14.0]

    if :recovered_unconf ∉ names(data)
        insertcols!(data, ncol(data) + 1, :recovered_unconf => zeros(Int, nrow(data)))
    end
    if :active_unconf ∉ names(data)
        insertcols!(data, ncol(data) + 1, :active_unconf => zeros(Int, nrow(data)))
    end

    model = SEIRModel(S, E, I, R, M, β, γ, α)
    model.data = data
    model
end

phuber_loss(a::Real, δ::Real) = δ^2 * (sqrt(1 + (a / δ)^2) - 1)
diff_phuber_loss(a::Real, δ::Real) = (sqrt(1 + (a / δ)^2)^-1) * a
phuber_loss(a::Vector{<:Real}, δ::Real) = δ^2 * (sqrt(1 + sum(a .* a) / δ^2) - 1)
diff_phuber_loss(a::Vector{<:Real}, δ::Real) = (sqrt(1 + (sum(a .* a) / δ^2))^-1) * a

function phuber_loss(δ::Real)
    function _phuber_loss(a)
        phuber_loss(a, δ)
    end
    _phuber_loss
end

function diff_phuber_loss(δ::Real)
    function _diff_phuber_loss(a)
        diff_phuber_loss(a, δ)
    end
    _phuber_loss
end

function estimate_unconf!(data::DataFrame, μ::Float64)
    for row ∈ eachrow(data)
        row.recovered_unconf = round((row.deaths / μ) - row.deaths - row.recovered)
        rate_conf = row.recovered / row.recovered_unconf
        row.active_unconf = round(row.active / rate_conf)
    end
    data
end

function loss(
    model::SEIRModel;
    loss_func::Function = phuber_loss(1.0),
    diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = paramsof(model)
    lethal_groups = model.lethal_groups
    data = dataof(model)
    sdata = to_dataframe(model)
    M = params[:M]
    μ = sum(M[lethal_groups]) / sum(M)
    estimate_unconf!(data, μ)
    #=
        rec_tot = loss_func(data.recovered .- mean(data.recovered))
        dth_tot = loss_func(data.deaths .- mean(data.deaths))
        act_tot = loss_func(data.active .- mean(data.active))
        rec_unconf_tot = loss_func(data.recovered_unconf .- mean(data.recovered_unconf))
        act_unconf_tot = loss_func(data.active_unconf .- mean(data.active_unconf))
        max_loss = rec_tot + dth_tot + act_tot + rec_unconf_tot + act_unconf_tot
    =#
    rec_loss = Float64[]
    act_loss = Float64[]
    dth_loss = Float64[]
    rec_unconf_loss = Float64[]
    act_unconf_loss = Float64[]

    diff_rec_loss = Float64[]
    diff_act_loss = Float64[]
    diff_dth_loss = Float64[]
    diff_rec_unconf_loss = Float64[]
    diff_act_unconf_loss = Float64[]

    for (i, si) ∈ zip(1:nrow(data), 1:4:nrow(sdata))
        push!(rec_loss, data.recovered[i] - sdata.recovered[si])
        push!(act_loss, data.active[i] - sdata.active[si])
        push!(dth_loss, data.deaths[i] - sdata.deaths[si])
        push!(rec_unconf_loss, data.recovered_unconf - sdata.recovered_unconf)
        push!(act_unconf_loss, data.active_unconf - sdata.active_unconf)

        push!(diff_rec_loss, sdata.diff_recovered[si])
        push!(diff_act_loss, sdata.diff_active[si])
        push!(diff_dth_loss, sdata.diff_deaths[si])
        push!(diff_rec_unconf_loss, sdata.diff_recovered_unconf)
        push!(diff_act_unconf_loss, sdata.diff_active_unconf)
    end
    loss_value = (
        loss_func(rec_loss) +
        loss_func(act_loss) +
        loss_func(dth_loss) +
        loss(rec_unconf_loss) +
        loss(act_unconf_loss)
    ) # / max_loss

    diff_loss_value = (
        diff_loss_func(rec_loss) .* diff_rec_loss +
        diff_loss_func(act_loss) .* diff_act_loss +
        diff_loss_func(dth_loss) .* diff_dth_loss +
        diff_loss_func(rec_unconf_loss) .* diff_rec_unconf_loss +
        diff_loss_func(act_unconf_loss) .* diff_act_unconf_loss
    )

    (loss_value, diff_loss_value)
end

function optimize_params(model::SEIRModel)
    lastparams = nothing
    diff_value = nothing
    function _calc_loss(params)
        paramsof!(model, unpack_params(params))
        (value, diff_value) = loss(model)
        lastparams = params
        value
    end
    function diff(params)
        if last_params != params
            _calc_loss(params)
        end
        return diff_value
    end

    n = model.ngroups
    params = paramsof(model)
    lower = zeros(Float64, 4 * n)
    maxM = fill(sum(params[:M]), n)
    maxβ = fill(10.0, n)
    maxαγ = fill(1.0, n)
    upper = pack_params((; M = maxM, β = maxβ, γ = maxαγ, α = maxαγ))
    minbox = Fminbox(NelderMead())
    optimize(
        _calc_loss,
        diff,
        lower,
        upper,
        pack_params(params),
        minbox;
        inplace = false,
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
    vectors = zip(M, β, γ, α)
    result = Float64[]
    for vec ∈ vectors
        append!(result, vec)
    end
    result
end

function pack_params(tup::Union{Tuple, NamedTuple})
    pack_params(tup...)
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

function unpack_params(P::Vector{Float64})
    Plen = length(P)
    np = 4
    Plen = length(P)
    M = view(P, 1:np:Plen)
    β = view(P, 2:np:Plen)
    γ = view(P, 3:np:Plen)
    α = view(P, 4:np:Plen)
    (; M = M, β = β, γ = γ, α = α)
end

function SEIR_ODE_fun(X, P, t)
    (S, E, I, R) = unpack_vars(X)
    (M, β, γ, α) = unpack_params(P)

    (dS, dE, dI, dR) = SEIRDeriv(SEIRVars(S, E, I, R), SEIRParams(M, β, γ, α))
    pack_vars(dS, dE, dI, dR)
end

function SEIR_ODEProblem(inivars::SEIRVars, params::SEIRParams)
    ODEProblem(SEIR_ODE_fun, pack_vars(inivars), (0.0, 180.0), pack_params(params))
end

function model_problem(model::SEIRModel)
    SEIR_ODEProblem(varsof(model), paramsof(model))
end

function model_solution(model::SEIRModel)
    solve(model_problem(model), Euler(); dt = 0.25)
end

function to_dataframe(problem::ODEProblem)
    solution = solve(problem, Euler(); dt = 0.25)
    days = solution.t
    n = length(days)
    ng = model.ngroups
    unconf = model.unconfirmed_groups
    lethal = model.lethal_groups
    exposed = zeros(Float64, n)
    active = zeros(Float64, n)
    active_unconf = zeros(Float64, n)
    recovered = zeros(Float64, n)
    recovered_unconf = zeros(Float64, n)
    deaths = zeros(Float64, n)
    diff_exposed = zeros(Float64, n)
    diff_active = zeros(Float64, n)
    diff_active_unconf = zeros(Float64, n)
    diff_recovered = zeros(Float64, n)
    diff_recovered_unconf = zeros(Float64, n)
    for i ∈ 1:n
        (S, E, I, R) = unpack_vars(solution[i])
        (dS, dE, dI, dR) = SEIRDeriv(S, E, I, R)
        exposed[i] = sum(E)
        active_unconf[i] = sum(I[unconf])
        active[i] = sum(I) - active_unconf[i]
        recovered_unconf[i] = sum(R[setdiff(unconf, lethal)])
        recovered[i] = sum(R[setdiff(1:ng, lethal)]) - recovered_unconf[i]
        deaths[i] = sum(R[lethal])
        diff_exposed[i] = sum(dE)
        diff_active_unconf[i] = sum(dI[unconf])
        diff_active[i] = sum(dI) - diff_active_unconf[i]
        diff_recovered_unconf[i] = sum(dR[setdiff(unconf, lethal)])
        diff_recovered[i] = sum(dR[setdiff(1:ng, lethal)]) - diff_recovered_unconf[i]
        diff_deaths[i] = sum(dR[lethal])
    end
    DataFrame(
        day = days,
        exposed = exposed,
        active = active,
        active_unconf = active_unconf,
        recovered = recovered,
        recovered_unconf = recovered_unconf,
        deaths = deaths,
        diff_exposed = exposed,
        diff_active_unconf = diff_active_unconf,
        diff_active = active,
        diff_recovered = recovered,
        diff_recovered_unconf = diff_recovered_unconf,
        diff_deaths = deaths,
    )
end

function to_dataframe!(model::SEIRModel)
    model.modeldata = to_dataframe(model_problem(model))
end

function to_dataframe(model::SEIRModel)
    to_dataframe(model_problem(model))
end

function model_plot(df::DataFrame, colnames = names(df); title = summary(df))
    xlabel = colnames[1]
    y1 = colnames[2]
    win = plot(
        df[!, xlabel],
        df[!, y1],
        xlabel = string(xlabel),
        label = string(y1),
        title = title,
        w = 3,
        legend = :right,
    )
    for yn ∈ colnames[3:end]
        plot!(
            df[!, xlabel],
            df[!, yn];
            xlabel = string(xlabel),
            label = string(yn),
            w = 3,
        )
    end
    win
end

function model_plot(model::SEIRModel)
    df = to_dataframe!(model)
    model_plot(df)
end

function model_plot(problem::ODEProblem)
    df = to_dataframe(problem)
    model_plot(df)
end
