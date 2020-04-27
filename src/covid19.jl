
@deftable Brazil(per_country_group) per_country_group[:Brazil]

@deftable per_state2(per_state, Brazil) begin
    nrec = Union{Int, Missing}[]
    for row ∈ eachrow(per_state)
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
    insertcols!(per_state, 5, :recovered => nrec)
end

@deftable per_city2(per_city, Brazil) begin
    nrec = Union{Int, Missing}[]
    for row ∈ eachrow(per_city)
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
    insertcols!(per_city, 6, :recovered => nrec)
end

@defcolumn closed(recovered, deaths) recovered .+ deaths
@defcolumn active(confirmed, closed) confirmed .- closed
@defgroup world2(world) (country, date) [
    (latitude, longitude) => mean,
    (confirmed, recovered, deaths, closed, active) => sum,
]

@deftable world3(world2, world_population, total_tests) begin
    df = join(world2, world_population; on = :country, kind = :left)
    gdf = join(df, total_tests; on = [:country, :date], kind = :left)
    keys = [:country, :date, :latitude, :longitude]
    vals1 = [:estimated_population, :total_tests, :test_kind]
    vals2 = [:confirmed, :recovered, :deaths, :closed, :active]
    select(gdf, keys ∪ vals1 ∪ vals2)
end

@defgroup per_country(world3) (country,) [
    (date, latitude, longitude) => identity,
    (estimated_population, total_tests, test_kind) => identity,
    (confirmed, recovered, deaths, closed, active) => identity,
]

@defgroup total(per_country) (date,) [
    (latitude, longitude) => mean,
    (estimated_population, total_tests) => sum,
    (test_kind,) => maximum,
    (confirmed, recovered, deaths, closed, active) => sum,
]

@defcolumn μ_closed(deaths, closed) 1.0 .* deaths ./ closed
@defcolumn μ_confirmed(deaths, confirmed) 1.0 .* deaths ./ confirmed

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
    if csse != nothing && !haskey(data, :Brazil) && haskey(csse, :Brazil)
        br = csse[:Brazil]
        data[:Brazil] = br
        if brasil_io != nothing
            brasil_io[:Brazil] = br
        end
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
        column_closed,
        column_active,
        group_world2,
        copy_pop_tests_tables,
        table_world3,
        group_per_country,
        group_total,
        table_Brazil,
        copy_brazil_tables,
        table_per_state2,
        table_per_city2,
        column_μ_closed,
        column_μ_confirmed,
        column_confirmed_per_1mi,
        column_deaths_per_1mi,
        column_recovered_per_1mi,
        column_closed_per_1mi,
        column_active_per_1mi,
        column_tests_per_confirmed,
        column_tests_per_1mi,
    ]
    Database{Dict{Symbol, Any}}(sources, kwargs, funcs)
end
