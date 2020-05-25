# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

# Parameter limits

const MIN_μ = 0.001
const MAX_μ = 0.5
const MIN_α = 0.01
const MAX_α = 0.5
const MIN_γ = 0.01
const MAX_γ = 1.0
const MIN_β = 0.01
const MAX_β = 0.5

const SEIR_σ_NAMES = (:σ_β, :σ_γ, :σ_α, :σ_μ, :σ_E)
const SEIRσ = NamedTuple{SEIR_σ_NAMES, NTuple{5, Float}}

SEIRσ(σ_β, σ_γ, σ_α, σ_μ, σ_E) = SEIRσ((σ_β, σ_γ, σ_α, σ_μ, σ_E))

phuber_loss(a::Real, δ::Real) = δ^2 * (sqrt(1 + (a / δ)^2) - 1)
diff_phuber_loss(a::Real, δ::Real) = (sqrt(1 + (a / δ)^2)^-1) * a
phuber_loss(a::Vector{<:Real}, δ::Real) = δ^2 * (sqrt(1 + sum(a .* a) / δ^2) - 1)
diff_phuber_loss(a::Vector{<:Real}, δ::Real) = (sqrt(1 + (sum(a .* a) / δ^2))^-1) * a

l2_loss(a::Real) = a^2
diff_l2_loss(a::Real) = 2 * a
l2_loss(a::Vector{<:Real}) = sum(a .* a)
diff_l2_loss(a::Vector{<:Real}) = 2 * a

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

function model_loss(
    model::SEIRModel;
    loss_func::Function = phuber_loss(1.0),
    kwargs...,
    # diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = parameters(model)
    ng = model.ngroups
    data = realdata(model)
    initial_date = default_kwarg(model, :initial_date, today() - Day(15))
    ini = findfirst(isequal(initial_date), data.date)
    if isnothing(ini)
        throw(ErrorException(
            "Initial date not found in model data",
            model,
            initial_date,
        ))
    end
    iend = nrow(data)
    data = data[ini:iend, :]
    ndays = nrow(data)
    mdata = modeldata(model)
    nrow(mdata) < ndays && return Inf

    @debug(
        "Computing maximum loss for model",
        model,
        recovered = Tuple(data.recovered),
        deaths = Tuple(data.deaths),
        active = Tuple(data.active),
    )
    rec_tot = loss_func(data.recovered .- mean(data.recovered))
    dth_tot = loss_func(data.deaths .- mean(data.deaths))
    act_tot = loss_func(data.active .- mean(data.active))
    max_loss = rec_tot + dth_tot + act_tot

    rec_loss = Float[]
    act_loss = Float[]
    dth_loss = Float[]
    for i ∈ 1:ndays
        push!(rec_loss, data.recovered[i] - mdata.recovered[i])
        push!(act_loss, data.active[i] - mdata.active[i])
        push!(dth_loss, data.deaths[i] - mdata.deaths[i])
    end
    loss_value = loss_func(rec_loss) + loss_func(act_loss) + loss_func(dth_loss)
    loss_value / max_loss
end

function optimize_parameters!(
    model::SEIRModel;
    packed_params = (:β, :γ, :α, :μ),
    kwargs...,
)
    n = model.ngroups
    lastparams = nothing
    model_params = parameters(model)
    function _calc_loss(arg)
        params = arg[1:(end - n)]
        newE = arg[(end - n + 1):end]
        newparams = unpack_params(
            params,
            model.ngroups;
            packed_params = packed_params,
            default_params = model_params,
        )

        parameters!(model, newparams)
        (S, E, I, R) = variables(model)
        variables!(model, SEIRVariables(S, newE, I, R))
        model_loss(model; kwargs...)
    end
    #= function diff(params)
        if lastparams != params
            _calc_loss(params)
        end
        return diff_value
    end =#

    maxM = fill(1.01 * sum(model_params[:M]), n)
    maxp = (; M = maxM, β = MAX_β, γ = MAX_γ, α = MAX_α, μ = MAX_μ)
    minp = (; M = 0.1 * maxM, β = MIN_β, γ = MIN_γ, α = MIN_α, μ = MIN_μ)
    maxE, minE = maxM ./ 2, fill(0.0, n)
    upper = pack_params(maxp; packed_params = packed_params)
    append!(upper, maxE)
    lower = pack_params(minp; packed_params = packed_params)
    append!(lower, minE)
    minbox = Fminbox(NelderMead())
    E0 = variables(model)[:E]
    opt_params = pack_params(model_params; packed_params = packed_params)
    append!(opt_params, Float.(E0))
    optimize(_calc_loss, lower, upper, opt_params, minbox)
end

optimize_parameters(model::SEIRModel; kwargs...) =
    optimize_parameters!(copy(model); kwargs...)

function optimize_parameters!(data::AbstractDict; kwargs...)
    for (key, subdata) ∈ data
        @debug "Optimizing parameters." key _debuginfo(subdata)
        optimize_parameters!(subdata; kwargs...)
    end
    data
end

function optimize_parameters(data::AbstractDict; kwargs...)
    result = empty(data)
    for (key, subdata) ∈ data
        @debug "Optimizing parameters." key _debuginfo(subdata)
        result[key] = optimize_parameters(subdata; kwargs...)
    end
    result
end

optimize_parameters(data::AbstractDatabase; kwargs...) =
    optimize_parameters(modeldict(data); kwargs...)

function optimize_parameters!(data::AbstractDatabase; kwargs...)
    result = optimize_parameters!(modeldict(data); kwargs...)
    parameters!(data)
    result
end

optimize_parameters(data; kwargs...) = data
optimize_parameters!(data; kwargs...) = data

function search_parameters(model::AbstractEndemicModel) end

function estimate_μ(data::AbstractDataFrame; ndays = 14, kwargs...)
    deaths, recovered = data.deaths, data.recovered
    deaths[end] == 0 && return [NaN]
    ind1 = findfirst(.!ismissing.(deaths) .& (deaths .> 0))
    ind2 = findfirst(.!ismissing.(recovered) .& (recovered .>= 100))
    (isnothing(ind1) || isnothing(ind2)) && return [NaN]

    ind = max(ind1, ind2, length(deaths) - ndays)
    dth, rec = deaths[ind:end], recovered[ind:end]
    vals = filter(x -> !ismissing(x) && isfinite(x) && x > 0, dth ./ (dth .+ rec))
    @debug(
        "Estimating μ = deaths / (deaths + recovered)",
        deaths = Tuple(dth),
        recovered = Tuple(rec),
        values = Tuple(vals)
    )
    isempty(vals) && return (NaN, Inf)
    val = mean(vals)
    if val >= MAX_μ || val <= MIN_μ
        return (NaN, Inf)
    else
        return (val, StatsBase.std(vals; mean = val))
    end
end

estimate_μ(model::AbstractEndemicModel; kwargs...) =
    estimate_μ(realdata(model); kwargs...)

function estimate_α(
    data::AbstractDataFrame;
    ndays = 7,
    μ_pair = estimate_μ(data),
    kwargs...,
)
    (μ, σ_μ) = μ_pair
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff_closed) .& :active .> 0)
    dR = data.diff_closed
    I = data.active
    isempty(I) && return (NaN, Inf)
    @debug "Estimating α = dR / I." dR = Tuple(dR) I = Tuple(I)
    vals = dR ./ I
    val = mean(vals)
    if val >= MAX_α || val <= MIN_α
        return (NaN, Inf)
    else
        return (val, StatsBase.std(vals; mean = val))
    end
end

estimate_α(model::AbstractEndemicModel; kwargs...) =
    estimate_α(realdata(model); kwargs...)

function _γ_root(d1, d2, d3, I, α)
    a = -I .* d2 .+ d1 .^ 2 .- I .* α .* d1
    b = d1 .* d2 .- I .* α .* d2
    c = -I .* d3
    Δ = b .^ 2 - 4 .* a .* c

    vals = [
        ai < 0 ? (-sqrt(abs(Δi)) - bi) / (2ai) : (sqrt(abs(Δi)) - bi) / (2ai)
        for (ai, bi, Δi) ∈ zip(a, b, Δ)
    ]
    root = mean(vals)

    @debug(
        "Estimating γ as the mean of max roots of equations `a * γ^2 + b * γ + c = 0`.",
        a = Tuple(a),
        b = Tuple(b),
        c = Tuple(c),
        Δ = Tuple(Δ),
        values = Tuple(vals),
        root = root
    )

    (root, StatsBase.std(vals; mean = root))
end

function estimate_γ(
    data::AbstractDataFrame;
    ndays = 12,
    μ_pair = estimate_μ(data),
    α_pair = estimate_α(data; μ_pair = μ_pair),
    kwargs...,
)
    (μ, σ_μ) = μ_pair
    (α, σ_α) = α_pair
    # Try to find at least 4 valid entries
    ini = something(findlast(!ismissing, data.diff3_infected), 8) - 3
    # Get numbers from last `ndays` days
    ini = min(ini, max(1, nrow(data) - ndays))
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff3_infected) .& .!ismissing.(:active))

    d1 = data.diff_infected
    d2 = data.diff2_infected
    d3 = data.diff3_infected

    I = data.active
    isempty(I) && return (NaN, Inf)
    (val, σ_γ) = _γ_root(d1, d2, d3, I, α)
    if val >= MAX_γ
        return (0.99 * MAX_γ, σ_γ + val - MAX_γ)
    elseif val <= MIN_γ
        return (1.01 * MIN_γ, σ_γ + MIN_γ - val)
    else
        return (val, σ_γ)
    end
end

estimate_γ(model::AbstractEndemicModel; kwargs...) =
    estimate_γ(realdata(model); kwargs...)

function estimate_β(
    data::AbstractDataFrame;
    ndays = 14,
    μ_pair = estimate_μ(data),
    α_pair = estimate_α(data; μ_pair = μ_pair),
    γ_pair = estimate_γ(data; μ_pair = μ_pair, α_pair = α_pair),
    kwargs...,
)
    (γ, σ_γ) = γ_pair
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff2_infected))

    d1 = data.diff_infected
    d2 = data.diff2_infected
    I = data.active

    @debug(
        "Estimating β = (d²(I + R) + γ * d(I + R)) / (γ * I + d(I + R)).",
        d1 = Tuple(d1),
        d2 = Tuple(d2),
        I = Tuple(I),
        γ_pair = γ_pair,
    )
    vals = (d2 .+ γ .* d1) ./ (γ .* I .+ d1)
    isempty(vals) && return (NaN, Inf)
    val = mean(vals)
    σ_β = StatsBase.std(vals; mean = val)
    if val >= MAX_β
        return (0.99 * MAX_β, σ_β + val - MAX_β)
    elseif val <= MIN_β
        return (1.01 * MIN_β, σ_β + MIN_β - val)
    else
        return (val, σ_β)
    end

end

function estimate_exposed(
    data::AbstractDataFrame;
    μ_pair = estimate_μ(data; kwargs...),
    α_pair = estimate_α(data; μ_pair = μ_pair, kwargs...),
    γ_pair = estimate_γ(data; μ_pair = μ_pair, α_pair = α_pair, kwargs...),
    kwargs...,
)
    (μ, σ_μ) = μ_pair
    (α, σ_α) = α_pair
    (γ, σ_γ) = γ_pair
    d1 = data.diff_infected

    @debug "Estimating E = d(I + R) / γ." d1 = Tuple(d1) γ = Tuple(γ)
    E = d1 ./ γ
    E = [ismissing(elt) ? missing : isfinite(elt) ? round(Int, elt) : 0 for elt ∈ E]
    σ_E = σ_γ .* d1 ./ γ .^ 2
    (max.(E, 0), σ_E)
end

function estimate_exposed!(data::AbstractDataFrame; kwargs...)
    (E, σ_E) = estimate_exposed(data; kwargs...)
    if :exposed ∉ names(data)
        insertcols!(data, ncol(data) + 1, :exposed => E)
    else
        data.exposed = E
    end
    (E, σ_E)
end

estimate_exposed(model::AbstractEndemicModel; kwargs...) =
    estimate_exposed(realdata(model); kwargs...)
estimate_exposed!(model::AbstractEndemicModel; kwargs...) =
    estimate_exposed!(realdata(model); kwargs...)

function SEIRModel(
    data::AbstractDataFrame;
    minimum_infected_factor = option(:minimum_infected_factor),
    ndays = 14,
    kwargs...,
)
    @debug "Computing SEIR model for data" _debuginfo(data)
    numpeople = data.estimated_population[1]
    istart = findfirst(data.infected .>= numpeople * minimum_infected_factor)
    iend = findlast(!ismissing, data.recovered)
    data = data[istart:iend, :]

    idx = nrow(data) - ndays
    idx <= 0 && return missing

    initial_date = data.date[idx]
    M = Float[numpeople]
    μ_pair = estimate_μ(data; kwargs...)
    α_pair = estimate_α(data; μ_pair = μ_pair, kwargs...)
    γ_pair = estimate_γ(data; μ_pair = μ_pair, α_pair = α_pair, kwargs...)
    β_pair =
        estimate_β(data; μ_pair = μ_pair, α_pair = α_pair, γ_pair = γ_pair, kwargs...)

    (μ, σ_μ) = μ_pair
    (α, σ_α) = α_pair
    (γ, σ_γ) = γ_pair
    (β, σ_β) = β_pair

    (E_vec, σ_E_vec) = estimate_exposed!(
        data;
        μ_pair = μ_pair,
        α_pair = α_pair,
        γ_pair = γ_pair,
        kwargs...,
    )
    E = Float[E_vec[idx]]
    σ_E = Float(σ_E_vec[idx])
    I = Float[data.active[idx]]
    R = Float[data.closed[idx]]
    S = M - E - I - R
    vars, params = SEIRVariables(S, E, I, R), SEIRParameters(M, [β], [γ], [α], [μ])
    σ = SEIRσ(σ_β, σ_γ, σ_α, σ_μ, σ_E)

    SEIRModel(
        vars,
        params;
        realdata = data,
        initial_date = initial_date,
        σ_all = σ,
        kwargs...,
    )
end

_SEIRModel_exceptions = []

function SEIRModel!(destiny::AbstractDataDict, data::AbstractDict; kwargs...)
    function _try(data)
        try
            val = SEIRModel(data; kwargs...)
            val isa AbstractDataDict && isempty(val) ? missing : val
        catch exception
            exception isa InterruptException && rethrow(exception)
            @debug "Exception during computation of SEIR model." exception
            push!(_SEIRModel_exceptions, (data, exception))
            missing
        end
    end
    for (key, subdata) ∈ data
        val = _try(subdata)
        !ismissing(val) && (destiny[key] = val)
    end
    destiny
end

SEIRModel(data::AbstractDict; kwargs...) =
    SEIRModel!(empty(data, Symbol, Any), data; kwargs...)

SEIRModel(database::AbstractDatabase; kwargs...) =
    SEIRModel(datadict(database); kwargs...)

function SEIRModel!(database::AbstractDatabase; kwargs...)
    modeldict!(database, SEIRModel(database; kwargs...))
    parameters!(database; kwargs...)
end

function parameters!(destiny::AbstractDict, data::AbstractDict; kwargs...)
    columns = [param => OptFloat[] for param ∈ SEIR_PARAMS]
    σ_cols = [sym => OptFloat[] for sym ∈ SEIR_σ_NAMES]
    paramdf = DataFrame(
        :key => String[],
        :group_name => String[],
        columns...,
        σ_cols...,
        :R0 => OptFloat[],
        :loss_7_days => OptFloat[],
        :loss_14_days => OptFloat[],
    )
    # Auxiliar value
    infparams = SEIRσ(Inf, Inf, Inf, Inf, Inf)
    for (key, subdata) ∈ data
        @debug "Computing parameters." key _debuginfo(subdata)
        if subdata isa AbstractDict
            subdestiny = get!(destiny, key, empty(destiny))
            parameters!(subdestiny, subdata; kwargs...)
            continue
        end
        !(subdata isa AbstractEndemicModel) && continue
        n = ngroups(subdata)
        gnames = groupnames(subdata)
        param_pairs = pairs(parameters(subdata))
        loss7 = model_loss(subdata; ndays = 7, kwargs...)
        loss7 = isfinite(loss7) ? loss7 : missing
        loss14 = model_loss(subdata; ndays = 14, kwargs...)
        loss14 = isfinite(loss14) ? loss14 : missing

        for (i, gname) ∈ enumerate(gnames)
            σ_group = Symbol("σ_", gname)
            σ_pairs = pairs(default_kwarg(subdata, σ_group, infparams))
            params = [k => (isfinite(v[i]) ? v[i] : missing) for (k, v) ∈ param_pairs]
            σ_params = [k => (isfinite(v) ? v : missing) for (k, v) ∈ σ_pairs]
            R0 = param_pairs[:β][i] / param_pairs[:α][i]
            R0 = isfinite(R0) ? R0 : missing
            moreparams = (; R0 = R0, loss_7_days = loss7, loss_14_days = loss14)
            row = (;
                key = string(key),
                group_name = string(prettify(gname)),
                params...,
                σ_params...,
                moreparams...,
            )
            @debug "Computed params for SEIR model" row
            push!(paramdf, row)
        end
    end
    new_σ_cols = [Symbol("σ(", string(sym)[4:end], ")") for sym ∈ SEIR_σ_NAMES]
    rename!(paramdf, (SEIR_σ_NAMES .=> new_σ_cols)...)
    sort!(paramdf, :key)
    !isempty(paramdf) && (destiny[:MODEL_PARAMETERS] = paramdf)
    destiny
end

parameters(data::AbstractDict; kwargs...) =
    parameters!(empty(data, Symbol, Any), data; kwargs...)
parameters!(data::AbstractDataDict; kwargs...) = parameters!(data, data; kwargs...)

parameters(data::AbstractDatabase; kwargs...) = parameters(modeldict(data); kwargs...)
parameters!(data::AbstractDatabase; kwargs...) =
    parameters!(modeldict(data); kwargs...)
