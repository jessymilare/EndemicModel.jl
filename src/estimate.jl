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

_ratemean(itr) = geomean(Complex.(1.0 .+ itr)).re - 1.0

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
    ndays = nothing,
    kwargs...,
    # diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = parameters(model)
    ng = model.ngroups
    data = datadict(model)
    ndays = something(ndays, Dates.days(data.date[end] - data.date[1]))
    model_back = model_step(model, -ndays)
    istart, iend = nrow(data) - ndays, nrow(data)
    solution = model_solution(model; maxtime = ndays |> Float64)
    days = Int.(solution.t)
    n = length(days)
    M = params[:M]
    μ = params[:μ]

    rec_tot = loss_func(data.recovered .- mean(data.recovered))
    dth_tot = loss_func(data.deaths .- mean(data.deaths))
    act_tot = loss_func(data.active .- mean(data.active))
    max_loss = rec_tot + dth_tot + act_tot

    rec_loss = Float64[]
    act_loss = Float64[]
    dth_loss = Float64[]
    #=
    diff_rec_loss = Float64[]
    diff_act_loss = Float64[]
    diff_dth_loss = Float64[] =#
    for i ∈ 1:(iend - istart + 1)
        day = -Dates.days(data.date[end] - data.date[i])
        sol = solution[i]
        (S, E, I, R) = unpack_vars(sol)
        push!(rec_loss, data.recovered[i - 1 + istart] - sum((1 .- μ) .* R))
        push!(act_loss, data.active[i - 1 + istart] - sum(I))
        push!(dth_loss, data.deaths[i - 1 + istart] - sum(μ .* R))
        #=
        push!(diff_rec_loss, sdata.diff_recovered[si])
        push!(diff_act_loss, sdata.diff_active[si])
        push!(diff_dth_loss, sdata.diff_deaths[si]) =#
    end
    loss_value = loss_func(rec_loss) + loss_func(act_loss) + loss_func(dth_loss)
    #=
    diff_loss_value = (
        diff_loss_func(rec_loss) .* diff_rec_loss +
        diff_loss_func(act_loss) .* diff_act_loss +
        diff_loss_func(dth_loss) .* diff_dth_loss
    ) =#

    loss_value / max_loss
end

function optimize_params(
    model::SEIRModel;
    packed_params = (:β, :γ, :α, :μ),
    ndays = 14,
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
    append!(opt_params, Float64.(E0))
    optimize(_calc_loss, lower, upper, opt_params, minbox)
end

function estimate_μ(data::AbstractDataFrame; ndays = 14, kwargs...)
    [data.μ_closed_est[1]]
end

estimate_μ(model::AbstractEndemicModel; kwargs...) =
    estimate_μ(datadict(model); kwargs...)

function estimate_α(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    kwargs...,
)
    data = @where(data, .!ismissing.(:diff_closed))
    ini = max(1, nrow(data) - ndays)
    # Get numbers from last `ndays` days
    dR = data.diff_closed[ini:end]
    I = data.active[ini:end]
    [max(_ratemean(dR ./ I), 0.0)]
end

estimate_α(model::AbstractEndemicModel; kwargs...) =
    estimate_α(datadict(model); kwargs...)

const MIN_γ = 0.01

function _γ_root(d1, d2, d3, I, α)
    a = -I .* d2 .+ d1 .^ 2 .- I .* α .* d1
    b = d1 .* d2 .- I .* α .* d2
    c = -I .* d3
    Δ = b .^ 2 - 4 .* a .* c

    root1 = _ratemean((-sqrt.(abs.(Δ)) .- b) ./ (2a))
    root2 = _ratemean((sqrt.(abs.(Δ)) .- b) ./ (2a))

    max(root1, root2, 0.01)
end

function estimate_γ(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    kwargs...,
)
    data = @where(data, .!ismissing.(:diff3_confirmed) .& .!ismissing.(:active))
    # Get numbers from last `ndays` days
    ini = max(1, nrow(data) - ndays)
    d1 = data.diff_confirmed[ini:end]
    d2 = data.diff2_confirmed[ini:end]
    d3 = data.diff3_confirmed[ini:end]

    I = data.active[ini:end]

    [_γ_root(d1, d2, d3, I, α)]
end

estimate_γ(model::AbstractEndemicModel; kwargs...) =
    estimate_γ(datadict(model); kwargs...)

const MIN_β = 0.01

function estimate_β(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    γ = estimate_γ(data; μ = μ, α = α),
    kwargs...,
)
    data = @where(data, .!ismissing.(:diff2_confirmed))
    ini = max(1, nrow(data) - ndays)
    # Get numbers from last `ndays` days
    d1 = data.diff_confirmed[ini:end]
    d2 = data.diff2_confirmed[ini:end]

    I = data.active[ini:end]

    # `factor` is `MIN_γ` when `γ` is `MIN_γ`
    # and approximately 1.0 when `γ` is far from `MIN_γ`
    factor = sqrt.(γ .^ 2 .- MIN_γ^2) ./ γ .+ MIN_γ
    [max(_ratemean(factor .* (d2 .+ γ .* d1) ./ (γ .* I .+ d1)), MIN_β)]
end

function estimate_exposed(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data; kwargs...),
    α = estimate_α(data; μ = μ, kwargs...),
    γ = estimate_γ(data; μ = μ, α = α, kwargs...),
    kwargs...,
)
    d1 = data.diff_confirmed

    factor = sqrt.(γ .^ 2 .- MIN_γ^2) ./ γ .+ MIN_γ
    E = factor .* d1 ./ γ
    E = [ismissing(elt) ? missing : round(Int, elt) for elt ∈ E]
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
    estimate_exposed(datadict(model); kwargs...)
estimate_exposed!(model::AbstractEndemicModel; kwargs...) =
    estimate_exposed!(datadict(model); kwargs...)

function SEIRModel(
    data::AbstractDataFrame;
    minimum_confirmed_factor = option(:minimum_confirmed_factor),
    kwargs...,
)
    numpeople = data.estimated_population[1]
    data = data[data.confirmed .>= numpeople * minimum_confirmed_factor, :]
    M = Float64[numpeople]
    μ = estimate_μ(data; kwargs...)
    α = estimate_α(data; μ = μ, kwargs...)
    γ = estimate_γ(data; μ = μ, α = α, kwargs...)
    β = estimate_β(data; μ = μ, α = α, γ = γ, kwargs...)

    exposed = estimate_exposed!(data; μ = μ, α = α, γ = γ, kwargs...)
    idx = findlast(!ismissing, exposed)
    I = Float64[data.active[idx]]
    S = M - I
    R = Float64[data.closed[idx]]
    E = Float64[exposed[idx]]
    initial_date = data.date[idx]
    vars, params = SEIRVariables(S, E, I, R), SEIRParameters(M, β, γ, α, μ)

    SEIRModel(vars, params; data = data, initial_date = initial_date, kwargs...)
end

function SEIRModel(data::D; kwargs...) where {D <: AbstractDict}
    columns = [param => Vector{Float64}[] for param ∈ SEIR_PARAMS]
    paramdf = DataFrame(:key => String[], columns...)
    function _try(data)
        try
            model = SEIRModel(data; kwargs...)
            if model isa AbstractEndemicModel
                params = pairs(parameters(model))
                push!(paramdf, (; key = String(key), first.(params)...))
            end
        catch
            missing
        end
    end
    result = D(key => _try(subdata) for (key, subdata) ∈ data)
    !isempty(paramdf) && (result[:SEIR_model_parameters] = paramdf)
    filter!(p -> !ismissing(p[2]), result)
end

SEIRModel(database::AbstractDatabase; kwargs...) =
    SEIRModel(datadict(database); kwargs...)

function SEIRModel!(database::AbstractDatabase; kwargs...)
    model = modeldict!(database, SEIRModel(database; kwargs...))
    model
end
