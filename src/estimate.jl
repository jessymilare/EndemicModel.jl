# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

# Parameter limits

const MIN_μ = 0.001
const MAX_μ = 0.300

# Source: https://journals.lww.com/cmj/Fulltext/2020/05050/Persistence_and_clearance_of_viral_RNA_in_2019.6.aspx
# Average infeccious period of 9.5 (6.0 to 11.0) days
# Source: https://www.sciencedirect.com/science/article/pii/S0163445320301195
# Average infeccious period of 11.0 (10.0 to 12.0) days
const MIN_α = 1 / 12.0
const MAX_α = 1 / 6.0

# Source: https://www.acpjournals.org/doi/10.7326/M20-0504
# Incubation period 5.1 (4.5 to 5.8) days
const MIN_γ = 1 / 5.8
const MAX_γ = 1 / 4.5

const MIN_β = 0.001
const MAX_β = 0.5

const ε = 1e-2

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
    model_days = 15,
    kwargs...,
    # diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = parameters(model)
    ng = model.ngroups
    data = realdata(model)
    initial_date = default_kwarg(model, :initial_date, today() - Day(model_days))
    ini = findfirst(isequal(initial_date), data.date)
    if isnothing(ini)
        throw(ErrorException("Initial date not found in model data", model, initial_date))
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
    #rec_tot = loss_func(data.recovered .- mean(data.recovered))
    dth_tot = loss_func(data.deaths .- mean(data.deaths))
    act_tot = loss_func(data.active .- mean(data.active))
    max_loss = act_tot + dth_tot # + rec_tot

    #rec_loss = Float[]
    act_loss_vec = data.active .- mdata.active[1:ndays]
    dth_loss_vec = data.deaths .- mdata.deaths[1:ndays]
    #rec_loss_vec = data.recovered .- mdata.recovered[1:ndays]

    act_loss, dth_loss = loss_func(act_loss_vec), loss_func(dth_loss_vec)
    @debug(
        "Computing loss",
        Tuple(act_loss_vec),
        Tuple(dth_loss_vec),
        act_loss,
        dth_loss,
        max_loss
    )
    loss_value = act_loss + dth_loss
    loss_value / max_loss
end

function optimize_parameters!(
    model::SEIRModel,
    params = (:M, :β, :E);
    model_days = 30,
    kwargs...,
)
    packed_params = intersect(SEIR_PARAMETERS, params)
    n = model.ngroups
    lastparams = nothing
    model_params = parameters(model)
    model_vars = variables(model)

    function _calc_loss(arg)
        if :E ∈ params
            newparams = arg[1:(end-n)]
            newE = arg[(end-n+1):end]
        else
            newparams = arg
        end
        newparams = unpack_params(
            newparams,
            model.ngroups;
            packed_params = packed_params,
            default_params = model_params,
        )
        parameters!(model, newparams)
        if :E ∈ params || :M ∈ params
            (S, E, I, R) = variables(model)
            if :E ∈ params
                E = newE
            end
            M = newparams[:M]
            if :M ∈ params
                S = M .- E .- I .- R
            end
            variables!(model, SEIRVariables(S, E, I, R))
        end
        modeldata!(
            model,
            to_dataframe(model; diff_columns = false, maxtime = model_days + 1, kwargs...),
        )
        model_loss(model; kwargs...)
    end

    maxM = fill((1 + ε) * sum(model_params[:M]), n)
    minM = (1 + ε) * (model_vars[:I] .+ model_vars[:R])
    maxp = (; M = maxM, β = MAX_β, γ = MAX_γ, α = MAX_α, μ = MAX_μ)
    minp = (; M = minM, β = MIN_β, γ = MIN_γ, α = MIN_α, μ = MIN_μ)

    upper = pack_params(maxp; packed_params = packed_params)
    lower = pack_params(minp; packed_params = packed_params)
    opt_params = pack_params(model_params; packed_params = packed_params)

    if :E ∈ params
        maxE, minE = maxM ./ 2, fill(0.0, n)
        append!(upper, maxE)
        append!(lower, minE)
        E0 = model_vars[:E]
        append!(opt_params, Float.(E0))
    end
    @debug "Initial parameter values:" opt_params, lower, upper

    opt = optimize(
        _calc_loss,
        lower,
        upper,
        opt_params,
        Fminbox(NelderMead()),
        Optim.Options(x_tol = 1 / 256),
    )
    modeldata!(model, to_dataframe(model))
    opt
end

optimize_parameters(model::SEIRModel; kwargs...) =
    optimize_parameters!(copy(model); kwargs...)

function optimize_parameters!(data::AbstractDict; kwargs...)
    result = DataDict()
    for (key, subdata) ∈ data
        @debug "Optimizing parameters." key _debuginfo(subdata)
        result[key] = optimize_parameters!(subdata; kwargs...)
    end
    result
end

function optimize_parameters(data::AbstractDict; kwargs...)
    result = DataDict()
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

function optimize_variables!(model::SEIRModel; kwargs...)
    n = model.ngroups
    lastparams = nothing
    model_vars, model_params = variables(model), parameters(model)
    function _calc_loss(vars)
        newvars = unpack_vars(vars)

        variables!(model, newvars)
        model_loss(model; kwargs...)
    end

    maxS = (1 + ε) .* model_params[:M]
    maxv = (; S = maxS, E = maxS, I = maxS, R = maxS)
    minS = fill(-0.01, n)
    minv = (; S = minS, E = minS, I = minS, R = minS)
    upper = pack_vars(maxv)
    lower = pack_vars(minv)
    minbox = Fminbox(NelderMead())
    opt_vars = pack_vars(model_vars)
    optimize(_calc_loss, lower, upper, opt_vars, minbox)
end

optimize_variables(model::SEIRModel; kwargs...) =
    optimize_variables!(copy(model); kwargs...)

function search_parameters(model::AbstractEndemicModel) end

function estimate_μ(data::AbstractDataFrame; ndays = 7, kwargs...)
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
    if val >= 0.3 || val <= MIN_μ
        return NaN
    else
        return val
    end
end

estimate_μ(model::AbstractEndemicModel; kwargs...) = estimate_μ(realdata(model); kwargs...)

function estimate_α(data::AbstractDataFrame; ndays = 7, kwargs...)
    # Get numbers from last `ndays` days
    # ini = max(1, nrow(data) - ndays)
    # data = data[ini:end, :]
    # # Filter valid entries
    # data = @where(data, (.!ismissing.(:diff_closed)) .& (:active .> 0))
    # dR = data.diff_closed
    # I = data.active
    # isempty(I) && return (NaN, Inf)
    # @debug "Estimating α = dR / I." dR = Tuple(dR) I = Tuple(I)
    # vals = dR ./ I
    #
    # if any(val -> val >= 1 - ε, vals)
    #     return (MAX_α, StatsBase.std(vals; mean = MAX_α))
    # elseif any(val -> val <= ε, vals)
    #     return (MIN_α, StatsBase.std(vals; mean = MIN_α))
    # else
    #     val = min(MAX_α, max(MIN_α, geomean(vals .+ 1) - 1))
    #     return (val, StatsBase.std(vals; mean = val))
    # end

    # Source: https://journals.lww.com/cmj/Fulltext/2020/05050/Persistence_and_clearance_of_viral_RNA_in_2019.6.aspx
    # Average infeccious period of 9.5 (6.0 to 11.0) days
    # Source: https://www.sciencedirect.com/science/article/pii/S0163445320301195
    # Average infeccious period of 11.0 (10.0 to 12.0) days
    1 / 11.0
end

estimate_α(model::AbstractEndemicModel; kwargs...) = estimate_α(realdata(model); kwargs...)

# function _γ_root(d1, d2, d3, I, α)
#     a = -I .* d2 .+ d1 .^ 2 .- I .* α .* d1
#     b = d1 .* d2 .- I .* α .* d2
#     c = -I .* d3
#     Δ = b .^ 2 - 4 .* a .* c
#
#     vals = [
#         ai < 0 ? (-sqrt(abs(Δi)) - bi) / (2ai) : (sqrt(abs(Δi)) - bi) / (2ai)
#         for (ai, bi, Δi) ∈ zip(a, b, Δ)
#     ]
#
#     @debug(
#         "Estimating γ as the mean of max roots of equations `a * γ^2 + b * γ + c = 0`.",
#         a = Tuple(a),
#         b = Tuple(b),
#         c = Tuple(c),
#         Δ = Tuple(Δ),
#         values = Tuple(vals),
#     )
#
#     vals
# end

function estimate_γ(
    data::AbstractDataFrame;
    # ndays = 7,
    # α = estimate_α(data; kwargs...),
    kwargs...,
)
    # # Try to find at least 4 valid entries
    # ini = something(findlast(!ismissing, data.diff3_infected), 8) - 3
    # # Get numbers from last `ndays` days
    # ini = min(ini, max(1, nrow(data) - ndays))
    # data = data[ini:end, :]
    # # Filter valid entries
    # data = @where(data, .!ismissing.(:diff3_infected) .& .!ismissing.(:active))
    #
    # d1 = data.diff_infected
    # d2 = data.diff2_infected
    # d3 = data.diff3_infected
    #
    # I = data.active
    # isempty(I) && return (NaN, Inf)
    # vals = _γ_root(d1, d2, d3, I, α)
    #
    # if any(val -> val >= 1 - ε, vals)
    #     return (MAX_γ, StatsBase.std(vals; mean = MAX_γ))
    # elseif any(val -> val <= ε, vals)
    #     return (MIN_γ, StatsBase.std(vals; mean = MIN_γ))
    # else
    #     val = min(MAX_γ, max(MIN_γ, geomean(vals .+ 1) - 1))
    #     return (val, StatsBase.std(vals; mean = val))
    # end

    # Source: https://www.acpjournals.org/doi/10.7326/M20-0504
    # Incubation period 5.1 (4.5 to 5.8) days
    1 / 5.1
end

estimate_γ(model::AbstractEndemicModel; kwargs...) = estimate_γ(realdata(model); kwargs...)

function estimate_β(
    data::AbstractDataFrame;
    ndays = 7,
    α = estimate_α(data; kwargs...),
    γ = estimate_γ(data; α = α, kwargs...),
    kwargs...,
)
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff2_infected))

    if nrow(data) == 0
        return (mean([MIN_β, MAX_β]), (MAX_β - MIN_β) / 2)
    end

    d1 = data.diff_infected
    d2 = data.diff2_infected
    I = data.active
    M = data.estimated_population[1]
    R = data.closed

    @debug(
        "Estimating β = (d²(I + R) + γ * d(I + R)) / (γ * I + d(I + R)).",
        d1 = Tuple(d1),
        d2 = Tuple(d2),
        I = Tuple(I),
        γ = γ,
    )
    s_inv = M ./ (M .- I .- R)
    vals = s_inv .* (d2 .+ γ .* d1) ./ (γ .* I .+ d1)
    isempty(vals) && return (NaN, Inf)

    if any(val -> val >= 1 - ε, vals)
        return (1 - ε) * MAX_β
    elseif any(val -> val <= ε, vals)
        return (1 + ε) * MIN_β
    else
        val = (1 - ε) * MAX_β
        return val
    end

end

function estimate_exposed(
    data::AbstractDataFrame;
    μ = estimate_μ(data; kwargs...),
    α = estimate_α(data; kwargs...),
    γ = estimate_γ(data; α = α, kwargs...),
    kwargs...,
)
    d1 = data.diff_infected

    @debug "Estimating E = d(I + R) / γ." d1 = Tuple(d1) γ = Tuple(γ)
    E = d1 ./ γ
    E = [ismissing(elt) ? 0 : isfinite(elt) ? round(Int, elt) : 0 for elt ∈ E]
    max.(E, 0)
end

function estimate_exposed!(data::AbstractDataFrame; kwargs...)
    E = estimate_exposed(data; kwargs...)
    if "exposed" ∉ names(data)
        insertcols!(data, ncol(data) + 1, :exposed => E)
    else
        data.exposed = E
    end
    E
end

estimate_exposed(model::AbstractEndemicModel; kwargs...) =
    estimate_exposed(realdata(model); kwargs...)
estimate_exposed!(model::AbstractEndemicModel; kwargs...) =
    estimate_exposed!(realdata(model); kwargs...)

function SEIRModel(
    data::AbstractDataFrame;
    modeldata = nothing,
    minimum_infected = option(:minimum_infected),
    model_days = 30,
    optimize::Bool = true,
    kwargs...,
)
    if !isnothing(modeldata)
        return SEIRModel(
            data,
            modeldata;
            minimum_infected = minimum_infected,
            model_days = model_days,
            kwargs...,
        )
    end
    !hasproperty(data, :date) && return missing
    nrows = nrow(data)
    if any(data.date[i] + Day(1) != data.date[i+1] for i ∈ 1:nrows-1)
        # If dates are not sequential, give up
        return missing
    end

    @debug "Computing SEIR model for data" _debuginfo(data)

    numpeople = data.estimated_population[1]
    istart = findfirst(data.infected .>= minimum_infected)
    iend = findlast(!ismissing, data.recovered)
    (isnothing(istart) || isnothing(iend) || (iend <= istart)) && return missing
    data = data[istart:iend, :]

    idx = nrow(data) - model_days + 1
    idx <= 0 && return missing

    initial_date = data.date[idx]
    M = Float[numpeople]
    μ = estimate_μ(data; kwargs...)
    α = estimate_α(data; μ = μ, kwargs...)
    γ = estimate_γ(data; μ = μ, α = α, kwargs...)
    β = estimate_β(data; μ = μ, α = α, γ = γ, kwargs...)


    E_vec = estimate_exposed!(data; μ = μ, α = α, γ = γ, kwargs...)
    E = Float[E_vec[idx]]
    I = Float[data.active[idx]]
    R = Float[data.closed[idx]]
    S = M - E - I - R
    vars, params = SEIRVariables(S, E, I, R), SEIRParameters(M, [β], [γ], [α], [μ])

    model = SEIRModel(vars, params; realdata = data, initial_date = initial_date, kwargs...)

    if optimize
        @debug "Finding optimal parameters for model" model
        result = optimize_parameters!(model; model_days = model_days, kwargs...)
        @debug "Optimal parameters for model computed" result
    end

    model
end

function SEIRModel(
    data::AbstractDataFrame,
    modeldata::AbstractDataFrame;
    minimum_infected = option(:minimum_infected),
    kwargs...,
)
    initial_date = modeldata.date[1]
    M = data.estimated_population[1]

    E, I = modeldata.exposed, modeldata.active
    R = modeldata.recovered .+ modeldata.deaths
    S = M .- (E .+ I .+ R)

    dS = S[end] - S[end-1]
    dE = E[end] - E[end-1]
    dI = I[end] - I[end-1]
    dR = R[end] - R[end-1]

    β = mean(-dS / (I[end-1] * S[end-1] / M))
    γ = mean((dI + dR) / E[end-1])
    α = mean(dR / I[end-1])
    μ = modeldata.deaths[end] / (modeldata.recovered[end] + modeldata.deaths[end])

    SEIRModel(
        [S[1]],
        [E[1]],
        [I[1]],
        [R[1]],
        Float[M],
        [β],
        [γ],
        [α],
        [μ];
        initial_date = initial_date,
        realdata = data,
        modeldata = modeldata,
    )
end

_SEIRModel_exceptions = []

function SEIRModel!(
    destiny::AbstractDataDict,
    data::AbstractDict;
    modeldata = nothing,
    kwargs...,
)
    function _try(key, data)
        try
            submodel = isnothing(modeldata) ? nothing : get(modeldata, key, nothing)
            val = SEIRModel(data; modeldata = submodel, kwargs...)
            val isa AbstractDataDict && isempty(val) ? missing : val
        catch exception
            exception isa InterruptException && rethrow(exception)
            @debug "Exception during computation of SEIR model." exception
            push!(_SEIRModel_exceptions, (data, exception))
            missing
        end
    end
    for (key, subdata) ∈ data
        val = _try(key, subdata)
        !ismissing(val) && (destiny[key] = val)
    end
    destiny
end

SEIRModel(data::AbstractDict; kwargs...) =
    SEIRModel!(empty(data, Symbol, Any), data; kwargs...)

SEIRModel(database::AbstractDatabase; kwargs...) = SEIRModel(datadict(database); kwargs...)

function SEIRModel!(database::AbstractDatabase; kwargs...)
    modeldict!(database, SEIRModel(database; kwargs...))
    parameters!(database; kwargs...)
end

function parameters!(destiny::AbstractDict, data::AbstractDict; kwargs...)
    columns = [param => OptFloat[] for param ∈ SEIR_PARAMETERS]
    paramdf = DataFrame(
        :key => String[],
        :group_name => String[],
        columns...,
        :R0 => OptFloat[],
        :Rt => OptFloat[],
        :loss_7_days => OptFloat[],
        :loss_14_days => OptFloat[],
    )
    # Auxiliar value
    for (key, subdata) ∈ data
        @debug "Computing parameters." key _debuginfo(subdata)
        if subdata isa AbstractDict
            subdestiny = get!(destiny, key, empty(destiny))
            try
                parameters!(subdestiny, subdata; kwargs...)
            catch
            end
            continue
        end
        !(subdata isa AbstractEndemicModel) && continue

        n = ngroups(subdata)
        gnames = groupnames(subdata)
        params = pairs(parameters(subdata))
        vars = pairs(variables(subdata))
        loss7 = model_loss(subdata; ndays = 7, kwargs...)
        loss7 = isfinite(loss7) ? loss7 : missing
        loss14 = model_loss(subdata; ndays = 14, kwargs...)
        loss14 = isfinite(loss14) ? loss14 : missing

        E_last = subdata.realdata.exposed[end]
        I_last = subdata.realdata.infected[end]
        R_last = subdata.realdata.recovered[end]
        S_last = params[:M] .- E_last .- I_last .- R_last

        for (i, gname) ∈ enumerate(gnames)
            params = [k => (isfinite(v[i]) ? v[i] : missing) for (k, v) ∈ params]
            R0 = params[:β][i] / params[:α][i]
            R0 = isfinite(R0) ? R0 : missing
            Rt = ismissing(R0) ? R0 : R0 * S_last[i] / params[:M][i]
            moreparams = (; R0 = R0, Rt = Rt, loss_7_days = loss7, loss_14_days = loss14)
            row = (;
                key = string(key),
                group_name = string(prettify(gname)),
                params...,
                moreparams...,
            )
            @debug "Computed params for SEIR model" row
            push!(paramdf, row)
        end
    end
    sort!(paramdf, :key)
    !isempty(paramdf) && (destiny[:MODEL_PARAMETERS] = paramdf)
    destiny
end

parameters(data::AbstractDict; kwargs...) =
    parameters!(empty(data, Symbol, Any), data; kwargs...)
parameters!(data::AbstractDataDict; kwargs...) = parameters!(data, data; kwargs...)

parameters(data::AbstractDatabase; kwargs...) = parameters(modeldict(data); kwargs...)
parameters!(data::AbstractDatabase; kwargs...) = parameters!(modeldict(data); kwargs...)
