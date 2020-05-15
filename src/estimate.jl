# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

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
    ndays = 7,
    kwargs...,
    # diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = parameters(model)
    ng = model.ngroups
    data = realdata(model)
    ini = max(1, nrow(data) - ndays)
    ndays = nrow(data) - ini
    data = data[ini:end, :]
    model_back =
        model_step(model, -ndays; maxtime = Float(ndays), initial_date = data.date[1])
    mdata = modeldata(model_back)

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

function optimize_params(
    model::SEIRModel;
    packed_params = (:β, :γ, :α, :μ),
    ndays = 7,
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
        model_loss(model; ndays = ndays, kwargs...)
    end
    #= function diff(params)
        if lastparams != params
            _calc_loss(params)
        end
        return diff_value
    end =#

    maxM = fill(sum(model_params[:M]), n)
    maxαγ = fill(1.0, n)
    maxβ = fill(1.0, n)
    maxp = (; M = maxM, β = maxβ, γ = maxαγ, α = maxαγ, μ = maxαγ)
    maxE = maxM ./ 2
    upper = pack_params(maxp; packed_params = packed_params)
    append!(upper, maxE)
    lower = fill(-1.0e-16, length(upper))
    minbox = Fminbox(NelderMead())
    E0 = model.vars[:E]
    opt_params = pack_params(model_params; packed_params = packed_params)
    append!(opt_params, Float.(E0))
    optimize(_calc_loss, lower, upper, opt_params, minbox)
end

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
    isempty(vals) ? [NaN] : [mean(vals)]
end

estimate_μ(model::AbstractEndemicModel; kwargs...) =
    estimate_μ(realdata(model); kwargs...)

function estimate_α(
    data::AbstractDataFrame;
    ndays = 21,
    μ = estimate_μ(data),
    kwargs...,
)
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff_closed) .& :active .> 0)
    dR = data.diff_closed
    I = data.active
    isempty(I) && return [NaN]
    @debug "Estimating α = dR / I." dR = Tuple(dR) I = Tuple(I)
    [max(mean(dR ./ I), 0.0)]
end

estimate_α(model::AbstractEndemicModel; kwargs...) =
    estimate_α(realdata(model); kwargs...)

const MIN_γ = 0.01

function _γ_root(d1, d2, d3, I, α)
    a = -I .* d2 .+ d1 .^ 2 .- I .* α .* d1
    b = d1 .* d2 .- I .* α .* d2
    c = -I .* d3
    Δ = b .^ 2 - 4 .* a .* c

    root1 = mean((-sqrt.(abs.(Δ)) .- b) ./ (2a))
    root2 = mean((sqrt.(abs.(Δ)) .- b) ./ (2a))

    @debug(
        "Estimating γ as the mean of max roots of equations `a * γ^2 + b * γ + c = 0`.",
        a = Tuple(a),
        b = Tuple(b),
        c = Tuple(c),
        Δ = Tuple(Δ),
        root1,
        root2,
    )

    max(root1, root2)
end

function estimate_γ(
    data::AbstractDataFrame;
    ndays = 21,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    kwargs...,
)
    # Try to find at least 7 valid entries
    ini = something(findlast(!ismissing, data.diff3_confirmed) - 7, 1)
    # Get numbers from last `ndays` days
    ini = min(ini, max(1, nrow(data) - ndays))
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff3_confirmed) .& .!ismissing.(:active))

    d1 = data.diff_confirmed
    d2 = data.diff2_confirmed
    d3 = data.diff3_confirmed

    I = data.active
    isempty(I) && return [NaN]
    [min(max(_γ_root(d1, d2, d3, I, α), MIN_γ), 1.0)]
end

estimate_γ(model::AbstractEndemicModel; kwargs...) =
    estimate_γ(realdata(model); kwargs...)

const MIN_β = 0.01

function estimate_β(
    data::AbstractDataFrame;
    ndays = 14,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    γ = estimate_γ(data; μ = μ, α = α),
    kwargs...,
)
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    data = data[ini:end, :]
    # Filter valid entries
    data = @where(data, .!ismissing.(:diff2_confirmed))

    d1 = data.diff_confirmed
    d2 = data.diff2_confirmed
    I = data.active

    # `factor` is `MIN_γ` when `γ` is `MIN_γ`
    # and approximately 1.0 when `γ` is far from `MIN_γ`
    factor = sqrt.(γ .^ 2 .- MIN_γ^2) ./ γ .+ MIN_γ
    @debug(
        "Estimating β = (d²(I + R) + γ * d(I + R)) / (γ * I + d(I + R)).",
        d1 = Tuple(d1),
        d2 = Tuple(d2),
        I = Tuple(I),
        γ = Tuple(γ),
    )
    [max(factor .* (d2 .+ γ .* d1) ./ (γ .* I .+ d1) |> mean, MIN_β)]
end

function estimate_exposed(
    data::AbstractDataFrame;
    μ = estimate_μ(data; kwargs...),
    α = estimate_α(data; μ = μ, kwargs...),
    γ = estimate_γ(data; μ = μ, α = α, kwargs...),
    kwargs...,
)
    d1 = data.diff_confirmed

    factor = sqrt.(γ .^ 2 .- MIN_γ^2) ./ γ .+ MIN_γ
    @debug "Estimating E = d(I + R) / γ." d1 = Tuple(d1) γ = Tuple(γ)
    E = factor .* d1 ./ γ
    E = [ismissing(elt) ? missing : isfinite(elt) ? round(Int, elt) : 0 for elt ∈ E]
    max.(E, 0)
end

function estimate_exposed!(data::AbstractDataFrame; kwargs...)
    E = estimate_exposed(data; kwargs...)
    if :exposed ∉ names(data)
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
    minimum_confirmed_factor = option(:minimum_confirmed_factor),
    kwargs...,
)
    @debug "Computing SEIR model for data" _debuginfo(data)
    numpeople = data.estimated_population[1]
    data = data[data.confirmed .>= numpeople * minimum_confirmed_factor, :]
    M = Float[numpeople]
    μ = estimate_μ(data; kwargs...)
    α = estimate_α(data; μ = μ, kwargs...)
    γ = estimate_γ(data; μ = μ, α = α, kwargs...)
    β = estimate_β(data; μ = μ, α = α, γ = γ, kwargs...)

    exposed = estimate_exposed!(data; μ = μ, α = α, γ = γ, kwargs...)
    idx = findlast(!ismissing, exposed)
    I = Float[data.active[idx]]
    S = M - I
    R = Float[data.closed[idx]]
    E = Float[exposed[idx]]
    initial_date = data.date[idx]
    vars, params = SEIRVariables(S, E, I, R), SEIRParameters(M, β, γ, α, μ)

    SEIRModel(vars, params; realdata = data, initial_date = initial_date, kwargs...)
end

function SEIRModel!(destiny::AbstractDataDict, data::AbstractDict; kwargs...)
    columns = [param => OptFloat[] for param ∈ SEIR_PARAMS]
    paramdf = DataFrame(
        :key => String[],
        :group_name => String[],
        columns...,
        :R0 => OptFloat[],
    )
    function _try(key, data)
        try
            model = SEIRModel(data; kwargs...)
            model isa AbstractDict && return model

            gnames = string.(model.groupnames)
            param_pairs = pairs(parameters(model))
            for (i, gname) ∈ enumerate(gnames)
                params =
                    [k => (isfinite(v[i]) ? v[i] : missing) for (k, v) ∈ param_pairs]
                R0 = param_pairs[:β][i] / param_pairs[:α][i]
                @debug "Computed params for SEIR model" (; params...) R0
                push!(paramdf, (; key = key, group_name = gname, params..., R0 = R0))
            end
            model
        catch exception
            @debug "Exception during computation of SEIR model." exception
            missing
        end
    end
    for (key, subdata) ∈ data
        val = _try(string(key), subdata)
        !ismissing(val) && (destiny[key] = val)
    end
    !isempty(paramdf) && (destiny[:MODEL_PARAMETERS] = sort!(paramdf, :key))
    destiny
end

SEIRModel(data::AbstractDict; kwargs...) =
    SEIRModel!(empty(data, Symbol, Any), data; kwargs...)

SEIRModel(database::AbstractDatabase; kwargs...) =
    SEIRModel(datadict(database); kwargs...)

SEIRModel!(database::AbstractDatabase; kwargs...) =
    modeldict!(database, SEIRModel(database; kwargs...))
