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
const SEIR_PARAMS = (:M, :β, :γ, :α, :τ)
const SEIRVars = NamedTuple{SEIR_VARS, NTuple{4, Vector{Float64}}}
const SEIRDeriv = NamedTuple{SEIR_DERIV, NTuple{4, Vector{Float64}}}
const SEIRParams = NamedTuple{
    SEIR_PARAMS,
    Tuple{
        Vector{Float64},
        Array{Float64, 2},
        Vector{Float64},
        Vector{Float64},
        Array{Float64, 2},
    },
}

SEIRVars(S, E, I, R) = (SEIRVars((S, E, I, R)))
SEIRParams(M, β, γ, α, τ) = (SEIRParams((M, β, γ, α, τ)))

function SEIRDeriv(vars::SEIRVars, params::SEIRParams)
    (S, E, I, R) = vars
    (M, β, γ, α, τ) = params

    SEIRDeriv(S, E, I, R, M, β, γ, α, τ)
end

"""
Compute the derivative of a SEIR-type differential equation.
S, E, I, R must be vectors of the same length `n`.
M, γ must be vectors of length `n`.
β must be an matrix of length `n x n`.
τ must be an lower triangular matrix of length `n x n`, i.e., τ[i,j] = 0 for j >= i.

Equations:
M = S + E + I + R
dM = 0
dS = - β * (I ./ sum(M)) .* S
dE = β * (I ./ sum(M)) .* S - γ .* E
dI = γ .* E - α .* I + τ * I - τ' * I
dR = α .* I
"""
function SEIRDeriv(
    S::AbstractVector{Float64},
    E::AbstractVector{Float64},
    I::AbstractVector{Float64},
    R::AbstractVector{Float64},
    M::AbstractVector{Float64},
    β::AbstractArray{Float64, 2},
    γ::AbstractVector{Float64},
    α::AbstractVector{Float64},
    τ::AbstractArray{Float64, 2},
)
    StoE = β * (I ./ sum(M)) .* S
    EtoI = γ .* E
    ItoR = α .* I

    dS = max.(-StoE, -S)
    dE = max.(StoE - EtoI, -E)
    dI = max.(EtoI - ItoR + τ * I - τ' * I, -I)
    dR = max.(ItoR, -R)

    (; dS = dS, dE = dE, dI = dI, dR = dR)
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
    β::AbstractArray{Float64, 2},
    γ::AbstractVector{Float64},
    α::AbstractVector{Float64},
    τ::AbstractArray{Float64, 2},
    groupnames::Union{Nothing, Vector{Symbol}} = nothing,
)
    SEIRModel(
        (; S = S, E = E, I = I, R = R),
        (; M = M, β = β, γ = γ, α = α, τ = τ),
        groupnames,
    )
end

function SEIRModel(data::DataFrame)
    numpeople = data.estimated_population[1]
    deaths = data.deaths[end]
    recovered = data.recovered[end]
    μ = deaths / (recovered + deaths)
    M = numpeople * [1.0 - μ, μ, 0.0, 0.0]
    I = data.active[1] * [0.0, 0.0, 1.0 - μ, μ]
    S = M - I
    E = Float64[0.0, 0.0, 0.0, 0.0]
    R = Float64[0.0, 0.0, data.recovered[1], data.deaths[1]]
    β = [
        0.15 0.15 0.15 0.15
        0.15 0.15 0.15 0.15
        0.15 0.15 0.15 0.15
        0.15 0.15 0.15 0.15
    ]
    γ = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
    α = [1.0 / 17.0, 1.0 / 14.0, 1.0 / 17.0, 1.0 / 14.0]
    τ = [
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
        0.1 0.0 0.0 0.0 # 10% probability of confirmation - nonlethal cases
        0.0 0.9 0.0 0.0 # 90% probability of confirmation - lethal cases
    ]

    if :recovered_unconf ∉ names(data)
        insertcols!(data, ncol(data) + 1, :recovered_unconf => zeros(Int, nrow(data)))
    end
    if :active_unconf ∉ names(data)
        insertcols!(data, ncol(data) + 1, :active_unconf => zeros(Int, nrow(data)))
    end

    model = SEIRModel(S, E, I, R, M, β, γ, α, τ)
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
    _diff_phuber_loss
end

function estimate_unconf!(data::DataFrame, μ::Float64)
    μ <= 0 && return data
    for row ∈ eachrow(data)
        rec_unconf = round((row.deaths / μ) - row.deaths - row.recovered)
        row.recovered_unconf = rec_unconf
        rate_conf = rec_unconf == 0 ? 0 : row.recovered / rec_unconf
        row.active_unconf = rate_conf == 0.0 ? 0 : round(row.active / rate_conf)
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
    sdata = to_dataframe!(model)
    M = params[:M]
    μ = sum(M[lethal_groups]) / sum(M)
    estimate_unconf!(data, μ)

    rec_tot = loss_func(data.recovered .- mean(data.recovered))
    dth_tot = loss_func(data.deaths .- mean(data.deaths))
    act_tot = loss_func(data.active .- mean(data.active))
    rec_unconf_tot = loss_func(data.recovered_unconf .- mean(data.recovered_unconf))
    act_unconf_tot = loss_func(data.active_unconf .- mean(data.active_unconf))
    max_loss = rec_tot + dth_tot + act_tot + rec_unconf_tot + act_unconf_tot

    rec_loss = Float64[]
    act_loss = Float64[]
    dth_loss = Float64[]
    rec_unconf_loss = Float64[]
    act_unconf_loss = Float64[]
    #=
    diff_rec_loss = Float64[]
    diff_act_loss = Float64[]
    diff_dth_loss = Float64[]
    diff_rec_unconf_loss = Float64[]
    diff_act_unconf_loss = Float64[]
    =#
    if nrow(sdata) < 4 * nrow(data)
        return Inf
    end
    for i ∈ 1:nrow(data)
        si = 4 * (i - 1) + 1
        push!(rec_loss, data.recovered[i] - sdata.recovered[si])
        push!(act_loss, data.active[i] - sdata.active[si])
        push!(dth_loss, data.deaths[i] - sdata.deaths[si])
        push!(rec_unconf_loss, data.recovered_unconf[i] - sdata.recovered_unconf[si])
        push!(act_unconf_loss, data.active_unconf[i] - sdata.active_unconf[si])
        #=
        push!(diff_rec_loss, sdata.diff_recovered[si])
        push!(diff_act_loss, sdata.diff_active[si])
        push!(diff_dth_loss, sdata.diff_deaths[si])
        push!(diff_rec_unconf_loss, sdata.diff_recovered_unconf[si])
        push!(diff_act_unconf_loss, sdata.diff_active_unconf[si])
        =#
    end
    loss_value = (
        loss_func(rec_loss) +
        loss_func(act_loss) +
        loss_func(dth_loss) +
        loss_func(rec_unconf_loss) +
        loss_func(act_unconf_loss)
    )
    #=
    diff_loss_value = (
        diff_loss_func(rec_loss) .* diff_rec_loss +
        diff_loss_func(act_loss) .* diff_act_loss +
        diff_loss_func(dth_loss) .* diff_dth_loss +
        diff_loss_func(rec_unconf_loss) .* diff_rec_unconf_loss +
        diff_loss_func(act_unconf_loss) .* diff_act_unconf_loss
    )
    =#

    loss_value / max_loss
end

function optimize_params(model::SEIRModel, packed_params = (:γ, :α, :τ))
    lastparams = nothing
    model_params = paramsof(model)
    function _calc_loss(params)
        newparams = unpack_params(
            params,
            model.ngroups;
            packed_params = packed_params,
            default_params = model_params,
        )

        paramsof!(model, newparams)
        loss(model)
    end
    #=function diff(params)
        if lastparams != params
            _calc_loss(params)
        end
        return diff_value
    end=#

    n = model.ngroups
    params = paramsof(model)
    maxM = fill(sum(params[:M]), n)
    maxαγ = fill(1.0, n)
    maxβ = fill(1.0, (n, n))
    maxτ = fill(0.0, (n, n))
    for i ∈ 2:n
        for j ∈ 1:(i - 1)
            maxτ[i, j] = 1.0
        end
    end
    maxp = (; M = maxM, β = maxβ, γ = maxαγ, α = maxαγ, τ = maxτ)
    upper = pack_params(maxp; packed_params = packed_params)
    lower = fill(-0.001, length(upper))
    minbox = Fminbox(NelderMead())
    optimize(
        _calc_loss,
        #diff,
        lower,
        upper,
        pack_params(params; packed_params = packed_params),
        minbox,
        #inplace = false,
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

function pack_params(M, β, γ, α, τ)
    args = filter(!isnothing, (M, γ, α))
    vectors = zip(args...)
    result = Float64[]
    for v ∈ vectors
        append!(result, v)
    end
    if β !== nothing
        append!(result, vec(β))
    end
    if τ !== nothing
        n = size(τ, 1)
        for i ∈ 2:n
            for j ∈ 1:(i - 1)
                push!(result, τ[i, j])
            end
        end
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
        τ = :τ ∈ packed_params ? tup[:τ] : nothing
        pack_params(M, β, γ, α, τ)
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
    vecvars = [:M, :γ, :α] ∩ packed_params
    Plen = length(vecvars) * ngroups
    index = 1
    np = length(vecvars)
    if :M ∈ vecvars
        M = view(P, index:np:Plen)
        index += 1
    else
        M = default_params[:M]
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
    if :β ∈ packed_params
        βend = Plen + ngroups * ngroups
        β_vec = P[(Plen + 1):βend]
        β = Array{Float64, 2}(undef, ngroups, ngroups)
        for i ∈ 1:ngroups
            for j ∈ 1:ngroups
                β[i, j] = β_vec[(j - 1) * ngroups + i]
            end
        end
    else
        β = default_params[:β]
        βend = Plen
    end
    if :τ ∈ packed_params
        τ_vec = P[(βend + 1):end]
        τ = zeros(Float64, (ngroups, ngroups))
        k = 1
        for i ∈ 2:ngroups
            for j ∈ 1:(i - 1)
                τ[i, j] = τ_vec[k]
                k += 1
            end
        end
    else
        τ = default_params[:τ]
    end
    (; M = M, β = β, γ = γ, α = α, τ = τ)
end

function SEIR_ODE_fun(X, P, t)
    (S, E, I, R) = unpack_vars(X)
    n = length(S)
    (M, β, γ, α, τ) = unpack_params(P, n)

    (dS, dE, dI, dR) = SEIRDeriv(SEIRVars(S, E, I, R), SEIRParams(M, β, γ, α, τ))
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

function to_dataframe(problem::ODEProblem, model::SEIRModel)
    solution = solve(problem, Euler(); dt = 0.25)
    (M, β, γ, α, τ) = paramsof(model)
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
    diff_deaths = zeros(Float64, n)
    for i ∈ 1:n
        (S, E, I, R) = unpack_vars(solution[i])
        (dS, dE, dI, dR) = SEIRDeriv(S, E, I, R, M, β, γ, α, τ)
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
        diff_exposed = diff_exposed,
        diff_active_unconf = diff_active_unconf,
        diff_active = diff_active,
        diff_recovered = diff_recovered,
        diff_recovered_unconf = diff_recovered_unconf,
        diff_deaths = diff_deaths,
    )
end

function to_dataframe!(model::SEIRModel)
    model.modeldata = to_dataframe(model_problem(model), model)
end

function to_dataframe(model::SEIRModel)
    to_dataframe(model_problem(model), model)
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
        legend = :left,
    )
    for yn ∈ colnames[3:end]
        ynstr = string(yn)
        if length(ynstr) >= 5 && ynstr[1:5] != "diff_"
            plot!(
                win,
                df[!, xlabel],
                df[!, yn];
                xlabel = string(xlabel),
                label = ynstr,
                w = 3,
            )
        end
    end
    gui(win)
end

function model_plot(model::SEIRModel)
    df = to_dataframe!(model)
    model_plot(df)
end

function model_plot(problem::ODEProblem)
    df = to_dataframe(problem)
    model_plot(df)
end
