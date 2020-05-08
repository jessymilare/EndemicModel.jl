
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
    ndays = nothing,
    kwargs...,
    # diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = paramsof(model)
    ng = model.ngroups
    data = dataof(model)
    ndays = something(ndays, Dates.days(data.date[end] - data.date[1]))
    istart, iend = nrow(data) - ndays, nrow(data)
    solution = model_solution(model; maxtime = -ndays - 1 |> Float64)
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
        sol = solution[end - i + 1]
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
    model_params = paramsof(model)
    function _calc_loss(arg)
        params = arg[1:(end - n)]
        newE = arg[(end - n + 1):end]
        newparams = unpack_params(
            params,
            model.ngroups;
            packed_params = packed_params,
            default_params = model_params,
        )

        paramsof!(model, newparams)
        (S, E, I, R) = varsof(model)
        varsof!(model, SEIRVars(S, newE, I, R))
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

function estimate_μ(data::AbstractDataFrame; ndays = 14)
    [data.μ_closed_est[1]]
end

estimate_μ(model::AbstractModel; kwargs...) = estimate_μ(dataof(model); kwargs...)

function estimate_α(data::AbstractDataFrame; ndays = 7, μ = estimate_μ(data))
    dt = data[data.diff_closed .!== missing, :]
    # Get numbers from last `ndays` days
    dR = dt.diff_closed[(end - ndays):end]
    I = dt.active[(end - ndays):end]
    @show dR, I
    [geomean(dR ./ I)]
end

estimate_α(model::AbstractModel; kwargs...) = estimate_α(dataof(model); kwargs...)

function _γ_root(d1, d2, d3, I, α)
    Δ1 = 4 .* I .* d2 .^ 3
    Δ2 = (-3 .* d1 .^ 2 .+ 2 .* I .* α .* d1 .+ I .^ 2 .* α .^ 2) .* d2 .^ 2
    Δ3 = (2 .* I .^ 2 .* α .* d3 .- 6 .* I .* d3 .* d1) .* d2
    Δ4 = 4 .* d3 .* d1 .^ 3 .- 4 .* I .* α .* d3 .* d1 .^ 2
    Δ5 = I .^ 2 .* d3 .^ 2
    Δ = Δ1 + Δ2 + Δ3 + Δ4 + Δ5
    b = (I .* α .- d1) .* d2 .+ I .* d3
    a = I .* d2 .- d1 .^ 2 .+ I .* α .* d1

    root1 = geomean((-sqrt.(abs.(Δ)) .- b) ./ (2a))
    root2 = geomean((sqrt.(abs.(Δ)) .- b) ./ (2a))

    @show root1, root2
    max(root1, root2)
end

function estimate_γ(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
)
    dt = data[data.diff3_confirmed .!== missing, :]
    # Get numbers from last `ndays` days
    d1 = dt.diff_confirmed[(end - ndays):end]
    d2 = dt.diff2_confirmed[(end - ndays):end]
    d3 = dt.diff3_confirmed[(end - ndays):end]

    I = dt.active[(end - ndays):end]

    [_γ_root(d1, d2, d3, I, α)]
end

estimate_γ(model::AbstractModel; kwargs...) = estimate_γ(dataof(model); kwargs...)

function estimate_β(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    γ = estimate_γ(data; μ = μ, α = α),
)
    dt = data[data.diff3_confirmed .!== missing, :]
    # Get numbers from last `ndays` days
    d1 = dt.diff_confirmed[(end - ndays):end]
    d2 = dt.diff2_confirmed[(end - ndays):end]

    I = dt.active[(end - ndays):end]

    [geomean((d2 .+ γ .* d1) ./ (γ .* I .+ d1))]
end

function estimate_exposed(
    data::AbstractDataFrame;
    ndays = 7,
    μ = estimate_μ(data),
    α = estimate_α(data; μ = μ),
    γ = estimate_γ(data; μ = μ, α = α),
)
    d1 = data.diff_confirmed
    E = d1 ./ γ
    [ismissing(elt) ? missing : round(Int, elt) for elt ∈ E]
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

estimate_exposed(model::AbstractModel; kwargs...) =
    estimate_exposed(dataof(model); kwargs...)
estimate_exposed!(model::AbstractModel; kwargs...) =
    estimate_exposed!(dataof(model); kwargs...)

function SEIRModel(
    data::AbstractDataFrame;
    minimum_confirmed_factor = option(:minimum_confirmed_factor),
)
    numpeople = data.estimated_population[1]
    data = data[data.confirmed .>= numpeople * minimum_confirmed_factor, :]
    M = Float64[numpeople]
    μ = estimate_μ(data)
    α = estimate_α(data; μ = μ)
    γ = estimate_γ(data; μ = μ, α = α)
    β = estimate_β(data; μ = μ, α = α, γ = γ)

    estimate_exposed!(data; μ = μ, α = α, γ = γ)
    idx = findlast(!ismissing, data.exposed)
    I = Float64[data.active[idx]]
    S = M - I
    R = Float64[data.closed[idx]]
    E = Float64[data.exposed[idx]]

    model = SEIRModel(S, E, I, R, M, β, γ, α, μ; initial_date = data.date[idx])
    dataof!(model, data)
    model
end
