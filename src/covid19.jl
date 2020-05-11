# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

@defgroup world2(world) (country, date) [
    (latitude, longitude) => mean,
    (confirmed, recovered, deaths) => sum,
]

@deftable world3(world2, world_population, total_tests) begin
    df = join(world2, world_population; on = :country, kind = :left)
    gdf = join(df, total_tests; on = [:country, :date], kind = :left)
    keys = [:country, :date, :latitude, :longitude]
    vals1 = [:estimated_population, :total_tests, :test_kind]
    vals2 = [:confirmed, :recovered, :deaths]
    testkinds = filter(!ismissing, Set(gdf.test_kind))
    if length(testkinds) > 1
        testkind = if "samples tested" ∈ testkinds
            "samples tested"
        elseif "people tested" ∈ testkinds
            "people tested"
        else
            pop!(testkinds)
        end
        rows = map(elt -> ismissing(elt) || elt == testkind, gdf.test_kind)
        gdf = gdf[collect(rows), :]
    end
    sort!(gdf, [:country, :date])
    select(gdf, keys ∪ vals1 ∪ vals2)
end

@defgroup per_country(world3) (country,) [
    (date, latitude, longitude) => identity,
    (estimated_population, total_tests, test_kind) => identity,
    (confirmed, recovered, deaths) => identity,
]

@defgroup total(per_country) (date,) [
    (latitude, longitude) => mean,
    (estimated_population, total_tests) => sum,
    (test_kind,) => maximum,
    (confirmed, recovered, deaths) => sum,
]

@deftable states2(states, Brazil) begin
    nrec = Union{Int, Missing}[]
    for row ∈ eachrow(states)
        brrow = Brazil[Brazil.date .== row.date, :]
        push!(
            nrec,
            if !isempty(brrow)
                brrec = brrow.recovered[1]
                brconf = brrow.confirmed[1]
                conf = row.confirmed
                brconf == 0 ? 0 : round(Union{Int, Missing}, brrec * conf / brconf)
            else
                missing
            end,
        )
    end
    insertcols!(states, 5, :recovered => nrec)
end

@deftable cities2(cities, Brazil) begin
    nrec = Union{Int, Missing}[]
    for row ∈ eachrow(cities)
        brrow = Brazil[Brazil.date .== row.date, :]
        push!(
            nrec,
            if !isempty(brrow)
                brrec = brrow.recovered[1]
                brconf = brrow.confirmed[1]
                conf = row.confirmed
                brconf == 0 ? 0 : round(Union{Int, Missing}, brrec * conf / brconf)
            else
                missing
            end,
        )
    end
    insertcols!(cities, 6, :recovered => nrec)
end

@defgroup per_state(states2) (state,) [
    (date, city_ibge_code) => identity,
    (estimated_population,) => identity,
    (confirmed, recovered, deaths) => identity,
]

@defgroup per_city(cities2) (city,) [
    (state, date, city_ibge_code) => identity,
    (estimated_population,) => identity,
    (confirmed, recovered, deaths) => identity,
]

@defcolumn closed(recovered, deaths) recovered .+ deaths
@defcolumn active(confirmed, closed) confirmed .- closed

function _column_diff(column, date, ndays)
    n = length(date)
    result = Vector{Union{Float64, Missing}}(missing, n)
    half2 = floor(Int, ndays / 2)
    half1 = ndays - half2
    istart = 1 + half1
    iend = n - half2
    for i ∈ istart:iend
        if ( # consistency check
            date[i - half1] + Day(ndays) == date[i] + Day(half2) == date[i + half2] &&
            !ismissing(column[i - half1]) &&
            isfinite(column[i - half1]) &&
            column[i - half1] > 0
        )
            # Geometric mean of column[i+1]/column[i] over `ndays`
            diff_rate = abs(column[i + half2] / column[i - half1])^(1 / ndays) - 1
            if !ismissing(diff_rate) && isfinite(diff_rate)
                result[i] = diff_rate * column[i]
            end
        end
    end
    result
end

@defcolumn diff_confirmed(date, confirmed) _column_diff(confirmed, date, 7)
#= begin
    n = length(date)
    result = Vector{Union{Float64, Missing}}(missing, n)
    for i ∈ 5:(n - 3)
        if (
            date[i - 4] + Day(7) == date[i] + Day(3) == date[i + 3] && # consistency check
            !ismissing(confirmed[i - 4]) &&
            isfinite(confirmed[i - 4]) &&
            confirmed[i - 4] > 0
        )
            # Geometric mean of confirmed[i+1]/confirmed[i] over 7 days
            diff_rate = (confirmed[i + 3] / confirmed[i - 4])^(1 / 7)
            if !ismissing(diff_rate) && isfinite(diff_rate)
                result[i] = diff_rate * confirmed[i]
            end
        end
    end
    result
end
=#

@defcolumn diff_recovered(date, recovered) _column_diff(recovered, date, 7)
@defcolumn diff_deaths(date, deaths) _column_diff(deaths, date, 7)

@defcolumn diff_closed(diff_recovered, diff_deaths) diff_recovered .+ diff_deaths
@defcolumn diff_active(diff_confirmed, diff_closed) diff_confirmed .- diff_closed

@defcolumn diff2_confirmed(date, diff_confirmed) _column_diff(diff_confirmed, date, 5)
@defcolumn diff2_recovered(date, diff_recovered) _column_diff(diff_recovered, date, 5)
@defcolumn diff2_deaths(date, diff_deaths) _column_diff(diff_deaths, date, 5)

@defcolumn diff2_closed(diff2_recovered, diff2_deaths) diff2_recovered .+ diff2_deaths
@defcolumn diff2_active(diff2_confirmed, diff2_closed) diff2_confirmed .- diff2_closed

@defcolumn diff3_confirmed(date, diff2_confirmed) _column_diff(
    diff2_confirmed,
    date,
    3,
)
#=
@defcolumn diff3_recovered(date, diff2_recovered) _column_diff(
    diff2_recovered,
    date,
    3,
)
@defcolumn diff3_deaths(date, diff2_deaths) _column_diff(diff2_deaths, date, 3)

@defcolumn diff3_closed(diff3_recovered, diff3_deaths) diff3_recovered .+ diff3_deaths
@defcolumn diff3_active(diff3_confirmed, diff3_closed) diff3_confirmed .- diff3_closed
=#
@defcolumn diff4_confirmed(date, diff3_confirmed) _column_diff(
    diff3_confirmed,
    date,
    3,
)

@defcolumn diff5_confirmed(date, diff4_confirmed) _column_diff(
    diff4_confirmed,
    date,
    3,
)

@defcolumn diff6_confirmed(date, diff5_confirmed) _column_diff(
    diff5_confirmed,
    date,
    3,
)

const COVID_19_μ = 0.034

@defcolumn μ_closed(deaths, closed) begin
    replace(1.0 .* deaths ./ closed, NaN => missing)
end
@defcolumn μ_confirmed(deaths, confirmed) begin
    replace(1.0 .* deaths ./ confirmed, NaN => missing)
end

@defcolumn μ_closed_est(deaths, recovered) begin
    ind = findfirst(map(
        p -> begin
            (d, r) = p
            !ismissing(d) && !ismissing(r) && d >= 1 && r >= 1000
        end,
        zip(deaths, recovered),
    ))
    if deaths[end] == 0 || ind == nothing
        COVID_19_μ
    else
        ind = max(1, length(deaths) - 14, ind)
        dth, rec = deaths[ind:end], recovered[ind:end]
        vals = filter(x -> !ismissing(x) && isfinite(x) && x > 0, dth ./ (dth .+ rec))
        isempty(vals) ? COVID_19_μ : mean(vals)
    end
end

@defcolumn confirmed_per_1mi(confirmed, estimated_population) begin
    1.0e6 .* confirmed ./ estimated_population
end
@defcolumn deaths_per_1mi(deaths, estimated_population) begin
    1.0e6 .* deaths ./ estimated_population
end
@defcolumn recovered_per_1mi(recovered, estimated_population) begin
    1.0e6 .* recovered ./ estimated_population
end
@defcolumn closed_per_1mi(closed, estimated_population) begin
    1.0e6 .* closed ./ estimated_population
end
@defcolumn active_per_1mi(active, estimated_population) begin
    1.0e6 .* active ./ estimated_population
end

@defcolumn tests_per_confirmed(total_tests, confirmed) 1.0 .* total_tests ./ confirmed
@defcolumn tests_per_1mi(total_tests, estimated_population) begin
    1.0e6 .* total_tests ./ estimated_population
end

function copy_pop_tests_tables(data::AbstractDict)
    pop = get(data, :world_population, nothing)
    tests = get(data, :total_tests, nothing)
    csse = get(data, :csse, nothing)
    brasil_io = get(data, :brasil_io, nothing)
    if brasil_io != nothing
        if pop != nothing && !haskey(brasil_io, :world_population)
            brasil_io[:world_population] = pop
            updated = true
        end
    end
    if csse != nothing
        updated = false
        if pop != nothing && !haskey(csse, :world_population)
            csse[:world_population] = pop
            updated = true
        end
        if tests != nothing && !haskey(csse, :total_tests)
            csse[:total_tests] = tests
            true
        else
            updated
        end
    else
        for (_, subdata) ∈ data
            if subdata isa AbstractDict && copy_pop_tests_tables(subdata)
                return true
            end
        end
        false
    end
end

function copy_brazil_tables(data::AbstractDict)
    csse = get(data, :csse, nothing)
    brasil_io = get(data, :brasil_io, nothing)
    if (
        csse != nothing &&
        brasil_io != nothing &&
        !haskey(brasil_io, :Brazil) &&
        haskey(csse, :per_country_group)
    )
        brasil_io[:Brazil] = csse[:per_country_group][:Brazil]
        true
    else
        for (_, subdata) ∈ data
            if subdata isa AbstractDict && copy_brazil_tables(subdata)
                return true
            end
        end
        false
    end
end

function covid_19_database(
    sources = [:csse, :brasil_io, :world_population, :total_tests];
    kwargs...,
)
    funcs = [
        # Tables
        group_world2,
        copy_pop_tests_tables,
        table_world3,
        group_per_country,
        group_total,
        copy_brazil_tables,
        table_states2,
        table_cities2,
        group_per_state,
        group_per_city,
        # Columns
        column_closed,
        column_active,
        column_diff_confirmed,
        column_diff_recovered,
        column_diff_deaths,
        column_diff_closed,
        column_diff_active,
        column_diff2_confirmed,
        column_diff2_recovered,
        column_diff2_deaths,
        column_diff2_closed,
        column_diff2_active,
        column_diff3_confirmed,
        column_diff4_confirmed,
        column_diff5_confirmed,
        column_μ_closed,
        column_μ_confirmed,
        column_μ_closed_est,
        column_confirmed_per_1mi,
        column_deaths_per_1mi,
        column_recovered_per_1mi,
        column_closed_per_1mi,
        column_active_per_1mi,
        column_tests_per_confirmed,
        column_tests_per_1mi,
    ]
    db = Database{Dict{Symbol, Any}}(sources, kwargs, funcs)
end

function covid_19(; kwargs...)
    db = covid_19_database(; kwargs...)
    SEIRModel!(db; kwargs...)
    export_data(db; kwargs...)
    db
end
