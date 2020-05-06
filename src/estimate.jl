
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
    diff_loss_func::Function = diff_phuber_loss(1.0),
)
    params = paramsof(model)
    lethal_groups = model.lethal_groups
    data = dataof(model)
    ndays = Dates.days(data.date[end] - data.date[1]) + 240
    sdata = to_dataframe!(model; maxtime = Float64(ndays))
    M = params[:M]
    μ = sum(M[lethal_groups]) / sum(M)

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
    for i ∈ 1:nrow(data)
        push!(rec_loss, data.recovered[i] - sdata.recovered[i])
        push!(act_loss, data.active[i] - sdata.active[i])
        push!(dth_loss, data.deaths[i] - sdata.deaths[i])
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

function optimize_params(model::SEIRModel, packed_params = (:β, :γ, :α))
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
        model_loss(model)
    end
    #= function diff(params)
        if lastparams != params
            _calc_loss(params)
        end
        return diff_value
    end =#

    maxM = fill(sum(model_params[:M]), n)
    maxαγ = fill(1.0, n)
    maxβ = 1.0
    maxp = (; M = maxM, β = maxβ, γ = maxαγ, α = maxαγ)
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

function estimate_μ(data)
    data.μ_closed_est[1]
end

function estimate_α(data; ndays = 7, μ = estimate_μ(data))
    dt = data[data.diff_reccovered .!= missing, :]
    # Get numbers from last `ndays` days
    diff_rec = dt.diff_recovered[(end - ndays):ndays]
    diff_dth = dt.diff_deaths[(end - ndays):ndays]
    dR = [diff_red, diff_dth]
    active = dt.active[(end - ndays):ndays]
    I = active .* [μ, 1.0 - μ]

    mean.(I ./ dR)
end

function estimate_γ(data; ndays = 7, α = estimate_α(data), μ = estimate_μ(data))
    dt = data[data.diff3_confirmed .!= missing, :]
    # Get numbers from last `ndays` days
    d1 = dt.diff_confirmed[end - ndays, ndays]
    d2 = dt.diff2_confirmed[end - ndays, ndays]
    d3 = dt.diff3_confirmed[end - ndays, ndays]

    active = dt.active[end - ndays, ndays]
    I = active * [μ, 1.0 - μ]

    (
        sqrt.(
            4 .* I .* d2^3 + (-3 * d1^2 .+ 2 .* I .* α .* d1 .+ I^2 .* α^2) * d2^2 .+
            (2 .* I^2 .* α .* d3 .- 6 .* I .* d3 .* d1) * d2 .+ 4 .* d3 .* d1^3 .-
            4 .* I .* α .* d3 .* d1^2 .+ I .^ 2 .* d3^2,
        ) .+ (d1 .- I .* α) .* d2 .- I .* d3
    ) ./ (2 .* I .* d2 .- 2 .* d1^2 .+ 2 .* I .* α .* d1)
end

function estimate_β(
    data;
    ndays = 7,
    α = estimate_α(data),
    μ = estimate_μ(data),
    γ = estimate_γ(data),
)
    dt = data[data.diff3_confirmed .!= missing, :]
    # Get numbers from last `ndays` days
    d1 = dt.diff_confirmed[end - ndays, ndays]
    d2 = dt.diff2_confirmed[end - ndays, ndays]

    active = dt.active[end - ndays, ndays]
    I = active * [μ, 1.0 - μ]

    (d2 .+ γ .* d1) ./ (γ .* I .+ d1)
end

function SEIRModel(data::DataFrame)
    data = data[data.confirmed .>= 1000, :]
    μ = estimate_μ(data)
    M = numpeople * [(1.0 - μ), μ]
    I = [(1.0 - μ), μ] * data.active[1]
    S = M - I
    E = I
    R = Float64[data.recovered[1], data.deaths[1]]
    α = estimate_α(data; μ = μ)
    γ = estimate_γ(data; μ = μ, α = α)
    β = estimate_β(data; μ = μ, α = α, γ = γ)

    model = SEIRModel(S, E, I, R, M, β, γ, α, [:non_lethal, :lethal])
    model.data = data
    model
end
