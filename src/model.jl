const PARAMS = (:M, :βIa, :βIs, :βIc, :γIa, :γIs, :γIca, :γIcs, :αRs, :αRc)
const NPARAMS = length(PARAMS)
const VARIABLES = (:S, :E, :Ia, :Is, :Ic, :R, :Rc)
const NVARIABLES = length(VARIABLES)
const ParamsNamedTuple{N} =
    NamedTuple{PARAMS, NTuple{NPARAMS, Vector{N}}} where {N <: Real}
const VariablesNamedTuple{N} =
    NamedTuple{VARIABLES, NTuple{NVARIABLES, Vector{N}}} where {N <: Real}

"""
Compute a SEIR-type differential equation. `n` is the number of groups to
consider. The other arguments are expected to be vectors of length `n`.

Constant:

Functions in the differential equations:
S = susceptible
E = exposed (non-infeccious)
Ia = assymptomatic infeccious
Is = symptomatic infeccious
Ic = tested and confirmed infeccious
R = recovered / deaths
Rc = recovered / deaths after confirmation

Parameters:
M = total susceptible population of each group
βIa = contamination rate from assymptomatic
βIs = contamination rate from symptomatic
βIc = contamination rate from confirmed
γIa = E -> Ia rate
γIs = Ia -> Is rate
γIca = Ia -> Ic rate (test and confirmation of assymptomatic)
γIcs = Is -> Ic rate (test and confirmation of symptomatic)
αRs = Is -> R rate (recovery/death of symptomatic)
αRc = Ic -> R rate (recovery/death of tested/confirmed)

"""
function model_ode(
    n::Integer,
    S::AbstractVector{N},
    E::AbstractVector{N},
    Ia::AbstractVector{N},
    Is::AbstractVector{N},
    Ic::AbstractVector{N},
    R::AbstractVector{N},
    Rc::AbstractVector{N},
    M::AbstractVector{N},
    βIa::AbstractVector{N},
    βIs::AbstractVector{N},
    βIc::AbstractVector{N},
    γIa::AbstractVector{N},
    γIs::AbstractVector{N},
    γIca::AbstractVector{N},
    γIcs::AbstractVector{N},
    αRs::AbstractVector{N},
    αRc::AbstractVector{N},
)::NTuple{NVARIABLES, Vector{N}} where {N <: Real}

    β = sum((βIa .* Ia + βIs .* Is + βIc .* Ic) ./ M)
    StoE = β * S
    EtoIa = γIa .* E
    IatoIs = γIs .* Ia
    IatoIc = γIca .* Ia
    IstoIc = γIcs .* Is
    IstoR = αRs .* Is
    IctoR = αRc .* Ic

    dS = -StoE
    dE = StoE - EtoIa
    dIa = EtoIa - IatoIc - IatoIs
    dIs = IatoIs - IstoIc - IstoR
    dIc = IatoIc + IstoIc - IctoR
    dR = IstoR
    dRc = IctoR

    (dS, dE, dIa, dIs, dIc, dR, dRc)
end

function model_ode_function(n::Integer, M::AbstractVector{Float64})
    nv = NVARIABLES
    np = NPARAMS
    Xlen = nv * n
    Plen = np * n

    function _model_function(
        X::Vector{Float64},
        P::Vector{Float64},
        t::Float64,
    )::Vector{Float64}
        S = view(X, 1:nv:Xlen)
        E = view(X, 2:nv:Xlen)
        Ia = view(X, 3:nv:Xlen)
        Is = view(X, 4:nv:Xlen)
        Ic = view(X, 5:nv:Xlen)
        R = view(X, 6:nv:Xlen)
        Rc = view(X, 7:nv:Xlen)

        M = view(P, 1:np:Plen)
        βIa = view(P, 2:np:Plen)
        βIs = view(P, 3:np:Plen)
        βIc = view(P, 4:np:Plen)
        γIa = view(P, 5:np:Plen)
        γIs = view(P, 6:np:Plen)
        γIca = view(P, 7:np:Plen)
        γIcs = view(P, 8:np:Plen)
        αRs = view(P, 9:np:Plen)
        αRc = view(P, 10:np:Plen)

        (dS, dE, dIa, dIs, dIc, dR, dRc) = model_ode(
            n,
            S,
            E,
            Ia,
            Is,
            Ic,
            R,
            Rc,
            M,
            βIa,
            βIs,
            βIc,
            γIa,
            γIs,
            γIca,
            γIcs,
            αRs,
            αRc,
        )

        pack_variables(dS, dE, dIa, dIs, dIc, dR, dRc)
    end

    function _model_function(X::Vector{Float64}, P::Function, t)
        _model_function(X, P(t), t)
    end

    _model_function
end

function pack_variables(
    S::AbstractVector{N},
    E::AbstractVector{N},
    Ia::AbstractVector{N},
    Is::AbstractVector{N},
    Ic::AbstractVector{N},
    R::AbstractVector{N},
    Rc::AbstractVector{N},
)::Vector{N} where {N <: Real}
    vectors = zip(S, E, Ia, Is, Ic, R, Rc)
    result = N[]
    for vec ∈ vectors
        append!(result, vec)
    end
    result
end

function unpack_variables(
    n::Integer,
    X::AbstractVector{N},
)::VariablesNamedTuple{Float64} where {N <: Real}
    nv = NVARIABLES
    Xlen = nv * n
    S = view(X, 1:nv:Xlen)
    E = view(X, 2:nv:Xlen)
    Ia = view(X, 3:nv:Xlen)
    Is = view(X, 4:nv:Xlen)
    Ic = view(X, 5:nv:Xlen)
    R = view(X, 6:nv:Xlen)
    Rc = view(X, 7:nv:Xlen)
    (; S = S, E = E, Ia = Ia, Is = Is, Ic = Ic, R = R, Rc = Rc)
end

function pack_parameters(
    M::AbstractVector{N},
    βIa::AbstractVector{N},
    βIs::AbstractVector{N},
    βIc::AbstractVector{N},
    γIa::AbstractVector{N},
    γIs::AbstractVector{N},
    γIca::AbstractVector{N},
    γIcs::AbstractVector{N},
    αRs::AbstractVector{N},
    αRc::AbstractVector{N},
)::Vector{N} where {N <: Real}
    vectors = zip(M, βIa, βIs, βIc, γIa, γIs, γIca, γIcs, αRs, αRc)
    result = N[]
    for vec ∈ vectors
        append!(result, vec)
    end
    result
end

function unpack_parameters(
    n::Integer,
    P::AbstractVector{N},
)::ParamsNamedTuple{N} where {N <: Real}
    nv = NVARIABLES
    np = NPARAMS
    Plen = nv * n
    M = view(P, 1:np:Plen)
    βIa = view(P, 2:np:Plen)
    βIs = view(P, 3:np:Plen)
    βIc = view(P, 4:np:Plen)
    γIa = view(P, 5:np:Plen)
    γIs = view(P, 6:np:Plen)
    γIca = view(P, 7:np:Plen)
    γIcs = view(P, 8:np:Plen)
    αRs = view(P, 9:np:Plen)
    αRc = view(P, 10:np:Plen)
    (;
        M = M,
        βIa = βIa,
        βIs = βIs,
        βIc = βIc,
        γIa = γIa,
        γIs = γIs,
        γIca = γIca,
        γIcs = γIcs,
        αRs = αRs,
        αRc = αRc,
    )
end

# Aproximate values for covid-19
function model_parameters(
    num_people::Float64 = 1e6;
    mortality::Float64 = 0.034,
    daily_growth::Float64 = 0.15,
    incubation_days::Float64 = 5.5,
    incubation_infeccious_days::Float64 = 3.0,
    symptomatic_days_to_recovery::Float64 = 17.0,
    symptomatic_days_to_death::Float64 = 14.0,
    prob_assymptomatic_testing::Float64 = 0.30,
    prob_symptomatic_testing::Float64 = 0.65,
    prob_severe_case_testing::Float64 = 0.95,
    confirmed_isolation_factor::Float64 = 0.0,
    testing_numdays::Float64 = 3.0,
    kwargs...,
)::ParamsNamedTuple{Float64}
    @argcheck 0.0 <= mortality <= 1.0
    #@argcheck 0.0 <= daily_growth <= 1.0
    @argcheck 0.0 <= confirmed_isolation_factor <= 1.0

    M = num_people * [1.0 - mortality, mortality]
    βIa = 0.55 * [daily_growth, daily_growth]
    βIs = [daily_growth, daily_growth]
    βIc = (1 - confirmed_isolation_factor) * βIs

    gia = 1 / (incubation_days - incubation_infeccious_days)
    gis = 1 / (incubation_infeccious_days)
    γIa = [gia, gia]
    γIs = [gis, gis]
    γIca = (prob_assymptomatic_testing / incubation_infeccious_days) * [1.0, 1.0]
    γIcs = [
        prob_symptomatic_testing / testing_numdays,
        prob_severe_case_testing / testing_numdays,
    ]
    αRs = [(1.0) / symptomatic_days_to_recovery, (1.0) / symptomatic_days_to_death]
    αRc = [
        1.0 / (symptomatic_days_to_recovery - testing_numdays),
        1.0 / (symptomatic_days_to_death - testing_numdays),
    ]
    (;
        M = M,
        βIa = βIa,
        βIs = βIs,
        βIc = βIc,
        γIa = γIa,
        γIs = γIs,
        γIca = γIca,
        γIcs = γIcs,
        αRs = αRs,
        αRc = αRc,
    )
end

function model_initial_values(
    num_people::Float64 = 1e6;
    mortality::Float64 = 0.034,
    initial_infected_incubation::Float64 = 0.0,
    initial_infeccious_assymptomatic::Float64 = 1.0,
    initial_infeccious_symptomatic::Float64 = 0.0,
    initial_confirmed::Float64 = 0.0,
    initial_recovered::Float64 = 0.0,
    initial_confirmed_recovered::Float64 = 0.0,
    initial_confirmed_deaths::Float64 = 0.0,
    initial_deaths::Float64 = 0.0,
    kwargs...,
)::VariablesNamedTuple{Float64}
    mort_dist = [1.0 - mortality, mortality]
    S = num_people * mort_dist
    E = initial_infected_incubation * mort_dist
    Ia = initial_infeccious_assymptomatic * mort_dist
    Is = initial_infeccious_symptomatic * mort_dist
    Ic = initial_confirmed * mort_dist
    R = [initial_recovered, initial_deaths]
    Rc = [initial_confirmed_recovered, initial_confirmed_deaths]
    (; S = S, E = E, Ia = Ia, Is = Is, Ic = Ic, R = R, Rc = Rc)
end

function model_problem(
    num_people::Float64 = 1e6,
    maxtime::Float64 = 180.0;
    params = nothing,
    params_function = nothing,
    kwargs...,
)
    inivals = model_initial_values(num_people; kwargs...)
    M = inivals[1] .|> Float64
    func = model_ode_function(2, M)
    if params != nothing
        params = pack_parameters(params...)
    elseif params_function == nothing
        params = pack_parameters(model_parameters(num_people; kwargs...)...)
    else
        params = params_function
    end
    inivals_packed = pack_variables(inivals...)

    ODEProblem(func, inivals_packed, (0.0, maxtime), params)
end

function model_solution(problem::ODEProblem)
    solution = solve(problem, Tsit5())
end

function model_decode_solution(solution::ODESolution)
    days = solution.t
    n = length(days)
    infected = zeros(Float64, n)
    active_confirmed = zeros(Float64, n)
    recovered_total = zeros(Float64, n)
    deaths_total = zeros(Float64, n)
    recovered_confirmed = zeros(Float64, n)
    deaths_confirmed = zeros(Float64, n)
    for i ∈ 1:n
        (S, E, Ia, Is, Ic, R, Rc) = unpack_variables(2, solution[i])
        infected[i] = sum(E) + sum(Ia) + sum(Is) + sum(Ic)
        active_confirmed[i] = sum(Ic)
        recovered_total[i] = R[1] + Rc[1]
        deaths_total[i] = R[2] + Rc[2]
        recovered_confirmed[i] = Rc[1]
        deaths_confirmed[i] = Rc[2]
    end
    DataFrame(
        day = days,
        infected = infected,
        active_confirmed = active_confirmed,
        recovered_total = recovered_total,
        deaths_total = deaths_total,
        recovered_confirmed = recovered_confirmed,
        deaths_confirmed = deaths_confirmed,
    )
end

function model_decode_solution(problem::ODEProblem)
    solution = model_solution(problem)
    model_decode_solution(solution)
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

function model_plot(solution::ODESolution)
    df = model_decode_solution(solution)
    model_plot(df)
end

function model_plot(problem::ODEProblem)
    solution = model_solution(problem)
    df = model_decode_solution(solution)
    model_plot(df)
end

huber_loss(a::Real, δ::Real) = abs(a) <= δ ? a^2 / 2 : δ * (abs(a) - δ / 2)
huber_loss(a::Vector{<:Real}, δ::Real) =
    sum(abs.(a)) <= δ ? sum(a .* a) / 2 : δ * (sum(abs.(a)) - δ / 2)

function huber_loss(δ::Real)
    function _huber_loss(a)
        huber_loss(a, δ)
    end
    _huber_loss
end

phuber_loss(a::Real, δ::Real) = δ^2 * (sqrt(1 + (a / δ)^2) - 1)
phuber_loss(a::Vector{<:Real}, δ::Real) = δ^2 * (sqrt(1 + sum(a .* a) / δ^2) - 1)

function phuber_loss(δ::Real)
    function _phuber_loss(a)
        phuber_loss(a, δ)
    end
    _phuber_loss
end

function model_loss_function(
    data::DataFrame;
    loss_func::Function = phuber_loss(1.0),
    recovered_weight::Float64 = 1.0,
    deaths_weight::Float64 = 1.0,
    active_weight::Float64 = 1.0,
)
    initial_date::Date = data[1, :date]
    final_date::Date = data[end, :date]
    tcol = Float64[Dates.days(t - initial_date) for t ∈ data.date]
    rec_conf = data.recovered .|> Float64
    dth_conf = data.deaths .|> Float64
    act_conf = data.active .|> Float64
    len = nrow(data)

    function _loss_function(solution::ODESolution)
        if any((s.retcode != :Success for s in sol))
            Inf
        else
            soldata = model_decode_solution(solution)
            sol_day = soldata.day
            Δt = mean(sol_day[2:end] .- sol_day[1:(end - 1)])

            sol_rec_conf_itp = soldata.recovered_confirmed |> _interpolate
            sol_dth_conf_itp = soldata.deaths_confirmed |> _interpolate
            sol_act_conf_itp = soldata.active_confirmed |> _interpolate

            sol_rec_conf = zeros(Float64, len)
            sol_dth_conf = zeros(Float64, len)
            sol_act_conf = zeros(Float64, len)

            sol_idx = 1
            tini = sol_day[sol_idx]
            tend = sol_day[sol_idx + 1]

            for t ∈ 1:len
                while t > tend
                    sol_idx += 1
                    tini = sol_day[sol_idx]
                    tend = sol_day[sol_idx + 1]
                end
                Δt = tend - tini
                # Explanation: find t_idx such that sol_day[t_idx] == t
                # (i.e. supposing sol_day is linearly interpolated)
                # We know that, if t == tini, then t_idx = sol_idx
                # We know that, if t == tend, then t_idx = sol_idx + 1
                # The following linear formula satisfies the two conditions above
                t_idx = sol_idx * (tend - t) / Δt + (sol_idx + 1) * (t - tini) / Δt
                sol_rec_conf[t] = sol_rec_conf_itp(t_idx) |> Float64
                sol_dth_conf[t] = sol_dth_conf_itp(t_idx) |> Float64
                sol_act_conf[t] = sol_act_conf_itp(t_idx) |> Float64
            end

            rec_loss = loss_func(sol_rec_conf - rec_conf) * recovered_weight
            dth_loss = loss_func(sol_dth_conf - dth_conf) * deaths_weight
            act_loss = loss_func(sol_act_conf - act_conf) * active_weight

            rec_loss + dth_loss + act_loss
        end
    end
    _loss_function
end

_interpolate(col) = interpolate(col, BSpline(Quadratic(Free(OnGrid()))))

const _METAPARAMETERS = Dict{DataFrame, Dict}()

function model_parameters(data::DataFrame)
    sort!(data, :date)
    if hasproperty(data, :estimated_population) &&
       (population = data.estimated_population[1]) !== missing
        num_people = population |> Float64
        max_people = population |> Float64
    else
        num_people = 10.0 * data.confirmed[end]
        max_people = Inf
    end

    initial_date = data.date[1]
    incubation_days::Float64 = 5.5
    incubation_infeccious_days::Float64 = 3.0

    actini = data.active[1]
    actmax = maximum(data.active)
    Δact = actmax - actini
    row_at_actmax = findfirst(isequal(actmax), data.active)
    date_at_actmax = data.date[row_at_actmax]
    day_at_actmax = Dates.days(date_at_actmax - initial_date)
    closed_at_actmax = data.closed[row_at_actmax]
    deaths_at_actmax = data.deaths[row_at_actmax]

    closed = data.closed[end]
    deaths = data.deaths[end]
    μ_closed = closed == 0 ? 0.034 : Float64(deaths / closed)

    new_closed = [0.0]
    append!(new_closed, data.closed[2:end] - data.closed[1:(end - 1)])
    nclmax = maximum(new_closed)
    row_nclmax = findfirst(isequal(nclmax), new_closed)
    closed_at_nclmax = data.closed[row_nclmax]
    row_previous_open =
        findfirst(row -> row.active + row.closed >= closed_at_nclmax, eachrow(data))
    symptomatic_days_to_close = row_nclmax - row_previous_open + 1

    deaths_at_nclmax = data.deaths[row_nclmax]
    row_previous_infected = findfirst(
        row -> row.active >= 0.9 * (deaths_at_nclmax - row.deaths) / μ_closed,
        eachrow(data),
    )
    if row_previous_infected == nothing
        symptomatic_days_to_death = 14.0
    else
        symptomatic_days_to_death = row_nclmax - row_previous_infected + 1
    end

    daily_growth = day_at_actmax == 0 ? 0.15 : Δact^(1.0 / day_at_actmax)

    _METAPARAMETERS[data] = Dict(;
        num_people = num_people,
        closed = closed,
        deaths = deaths,
        actini = actini,
        actmax = actmax,
        row_at_actmax = row_at_actmax,
        nclmax = nclmax,
        closed_at_nclmax = closed_at_nclmax,
        deaths_at_nclmax = deaths_at_nclmax,
        mortality = μ_closed,
        daily_growth = daily_growth,
        incubation_days = Float64(incubation_days),
        incubation_infeccious_days = Float64(incubation_infeccious_days),
        symptomatic_days_to_recovery = Float64(symptomatic_days_to_close),
        symptomatic_days_to_death = Float64(symptomatic_days_to_death),
        symptomatic_days_to_close = Float64(symptomatic_days_to_close),
        prob_assymptomatic_testing = 0.30,
        prob_symptomatic_testing = 0.65,
        prob_severe_case_testing = 0.95,
        confirmed_isolation_factor = 0.0,
        testing_numdays = 3.0,
    )

    model_parameters(
        num_people;
        mortality = μ_closed,
        daily_growth = daily_growth,
        incubation_days = Float64(incubation_days),
        incubation_infeccious_days = Float64(incubation_infeccious_days),
        symptomatic_days_to_recovery = Float64(symptomatic_days_to_close),
        symptomatic_days_to_death = Float64(symptomatic_days_to_death),
        prob_assymptomatic_testing = 0.30,
        prob_symptomatic_testing = 0.65,
        prob_severe_case_testing = 0.95,
        confirmed_isolation_factor = 0.0,
        testing_numdays = 3.0,
    )
end

function model_initial_values(data::DataFrame)
    if !haskey(_METAPARAMETERS, data)
        model_parameters(data)
    end
    mparams = _METAPARAMETERS[data]
    num_people = mparams[:num_people]
    mortality = mparams[:mortality]

    initial_confirmed = data.active[1]
    prob_conf_symp = mparams[:prob_symptomatic_testing]
    prob_conf_assymp = mparams[:prob_assymptomatic_testing]
    prob_conf_total = (1 - prob_conf_symp) + (1 - prob_conf_assymp)

    initial_unconfirmed = initial_confirmed * prob_conf_total
    initial_infeccious_symptomatic = initial_unconfirmed / prob_conf_symp
    initial_infeccious_assymptomatic = initial_unconfirmed / prob_conf_assymp
    initial_confirmed_recovered = data.recovered[1]
    initial_confirmed_deaths = data.deaths[1]

    model_initial_values(
        num_people;
        mortality = mortality,
        initial_confirmed = initial_confirmed,
        initial_confirmed_recovered = initial_confirmed_recovered,
        initial_infeccious_symptomatic = initial_infeccious_symptomatic,
        initial_confirmed_deaths = initial_confirmed_deaths,
    )
end
