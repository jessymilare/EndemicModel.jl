# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

@defgroup world2(world) (country, date) [
    (latitude, longitude) => mean,
    (confirmed, recovered, deaths) => sum,
]

@deftable tests_no_duplicate(total_tests) begin
    groups = groupby(total_tests, [:country, :date])
    data = combine(groups) do subdf
        testkind = _choose_test(filter(!ismissing, subdf.test_kind))
        idx = findfirst(isequal(testkind), subdf.test_kind)
        subdf[idx, :]
    end
    data = sort!(data, [:country, :date])
end

function _choose_test(tests)
    if isempty(tests)
        missing
    elseif "people tested" ∈ tests
        "people tested"
    elseif "samples tested" ∈ tests
        "samples tested"
    elseif "units unclear" ∈ tests
        "units unclear"
    else
        @warn "Unrecognized test kind(s)." tests
        collect(tests)[1]
    end
end

@deftable world3(world2, world_population, world_gdp_per_capita, tests_no_duplicate) begin
    data = leftjoin(world2, world_population; on = :country)
    data = leftjoin(data, world_gdp_per_capita; on = [:country, :country_code])
    data = leftjoin(data, tests_no_duplicate; on = [:country, :date])
    rename!(data, :recovered => :recovered_confirmed)

    # Estimate recovered
    # Source: https://journals.lww.com/cmj/Fulltext/2020/05050/Persistence_and_clearance_of_viral_RNA_in_2019.6.aspx
    # Average infeccious period of 9.5 (6.0 to 11.0) days
    # Source: https://www.sciencedirect.com/science/article/pii/S0163445320301195
    # Average infeccious period of 11.0 (10.0 to 12.0) days
    infected_vec, countries = data.confirmed, data.country
    active, closed, last_infected = 0, 0, 0
    closed_vec = Int[]
    last_country = "VOID"

    for (infected, country) ∈ zip(infected_vec, countries)
        if country == last_country
            new_closed, new_active = active / 11.0, infected - last_infected
            closed += new_closed
            active += new_active - new_closed
            last_infected = infected
            push!(closed_vec, round(Int, closed))
        else
            active, closed, last_infected = 0, 0, 0
            last_country = country
            push!(closed_vec, 0)
        end
    end
    insertcols!(data, :recovered => closed_vec .- data.deaths)

    cols = [
        :country,
        :region,
        :income_group,
        :date,
        :latitude,
        :longitude,
        :estimated_population,
        :gdp_per_capita,
        :total_tests,
        :test_kind,
        :confirmed,
        :recovered_confirmed,
        :recovered,
        :deaths,
    ]
    select!(data, cols)
end

@defgroup per_country_aux(world3) (country, region) [
    (income_group,) => identity,
    (date, latitude, longitude) => identity,
    (estimated_population, gdp_per_capita) => identity,
    (total_tests, test_kind) => identity,
    (confirmed, recovered, deaths) => identity,
]

@defgroup per_country(per_country_aux) (country,) [
    (region, income_group) => identity,
    (date, latitude, longitude) => identity,
    (estimated_population, gdp_per_capita) => identity,
    (total_tests, test_kind) => identity,
    (confirmed, recovered, deaths) => identity,
]

@defgroup per_region_aux(world3) (region, date) [
    (estimated_population,) => sum,
    (total_tests,) => _sum ∘ skipmissing,
    (confirmed, recovered, deaths) => sum,
]

@defgroup per_region(per_region_aux) (region,) [
    (date,) => identity,
    (estimated_population,) => sum,
    (total_tests,) => _sum ∘ skipmissing,
    (confirmed, recovered, deaths) => sum,
]

@defgroup per_income_group_aux(world3) (income_group, date) [
    (estimated_population,) => sum,
    (total_tests,) => _sum ∘ skipmissing,
    (confirmed, recovered, deaths) => sum,
]

@defgroup per_income_group(per_income_group_aux) (income_group,) [
    (date,) => identity,
    (estimated_population,) => sum,
    (total_tests,) => _sum ∘ skipmissing,
    (confirmed, recovered, deaths) => sum,
]

@defgroup per_date(per_country) (date,) [
    (latitude, longitude) => mean,
    (estimated_population,) => sum,
    (total_tests,) => _sum ∘ skipmissing,
    (confirmed, recovered, deaths, infected, closed, active) => sum,
    (diff_infected, diff_recovered, diff_deaths, diff_closed, diff_active) => sum,
    (diff2_infected, diff3_infected, diff4_infected, diff5_infected) => sum,
]

@deftable per_date2(per_date, world_population, world_gdp_per_capita) begin
    wpop, wgdp = world_population, world_gdp_per_capita
    pop = wpop[wpop.country.=="World", :estimated_population]
    gdp = wgdp[wgdp.country.=="World", :gdp_per_capita]
    cols = [
        :date,
        :total_tests,
        :confirmed,
        :recovered,
        :deaths,
        :infected,
        :closed,
        :active,
    ]
    data = select(per_date, cols)
    insertcols!(data, 2, :estimated_population => pop[1])
    insertcols!(data, 3, :gdp_per_capita => gdp[1])
    insertcols!(data, 5, :test_kind => "units unclear")
    data
end

@deftable brazil2(brazil, total, estimates) begin
    keys = [:date, :estimated_population, :gdp_per_capita, :total_tests, :test_kind]
    br = select(brazil, keys)
    avg_ipc = estimates[estimates.state.=="total", :infected_per_confirmed][1]

    # Estimate recovered
    # Source: https://journals.lww.com/cmj/Fulltext/2020/05050/Persistence_and_clearance_of_viral_RNA_in_2019.6.aspx
    # Average infeccious period of 9.5 (6.0 to 11.0) days
    # Source: https://www.sciencedirect.com/science/article/pii/S0163445320301195
    # Average infeccious period of 11.0 (10.0 to 12.0) days
    infected_vec = round.(Int, avg_ipc * brazil.confirmed)
    active, closed, last_infected = 0, 0, 0
    closed_vec = Int[]
    for infected ∈ infected_vec
        new_closed, new_active = active / 11.0, infected - last_infected
        closed += new_closed
        active += new_active - new_closed
        last_infected = infected
        push!(closed_vec, round(Int, closed))
    end

    insertcols!(br, 1, :country => "Brazil")
    insertcols!(br, :confirmed => brazil.confirmed)
    insertcols!(br, :recovered_confirmed => brazil.recovered)
    insertcols!(br, :recovered => closed_vec .- brazil.deaths)
    insertcols!(br, :deaths => brazil.deaths)
    insertcols!(br, :infected => infected_vec)
    sort!(br, :date)
end

@deftable states2(states, brazil2) begin
    nrec = Union{Int,Missing}[]
    for row ∈ eachrow(states)
        brrow = brazil2[brazil2.date.==row.date, :]
        push!(
            nrec,
            if !isempty(brrow)
                brrec = brrow.recovered[1]
                brinfec = brrow.infected[1]
                infec = row.infected
                brinfec == 0 ? 0 : round(Union{Int,Missing}, brrec * infec / brinfec)
            else
                missing
            end,
        )
    end
    insertcols!(states, 5, :recovered => nrec)
end

@deftable cities2(cities, brazil2) begin
    nrec = Union{Int,Missing}[]
    for row ∈ eachrow(cities)
        brrow = brazil2[brazil2.date.==row.date, :]
        push!(
            nrec,
            if !isempty(brrow)
                brrec = brrow.recovered[1]
                brinfec = brrow.infected[1]
                infec = row.infected
                brinfec == 0 ? 0 : round(Union{Int,Missing}, brrec * infec / brinfec)
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
    (confirmed, recovered, deaths, infected) => identity,
]

@defgroup per_city(cities2) (city, state) [
    (state, date, city_ibge_code) => identity,
    (estimated_population,) => identity,
    (confirmed, recovered, deaths, infected) => identity,
]

@defcolumn infected(confirmed) confirmed
@defcolumn closed(recovered, deaths) recovered .+ deaths
@defcolumn active(infected, closed) infected .- closed

function _column_diff(column, date, ndays)
    @argcheck ndays >= 1
    n = length(date)
    result = Vector{OptFloat}(missing, n)
    half2 = floor(Int, ndays / 2)
    half1 = ndays - half2
    istart = 1 + half1
    iend = n - half2
    for i ∈ istart:iend
        if ( # consistency check
            date[i-half1] + Day(ndays) == date[i] + Day(half2) == date[i+half2] &&
            !ismissing(column[i-half1]) &&
            isfinite(column[i-half1])
        )
            result[i] = (column[i+half2] - column[i-half1]) / ndays
        end
    end
    result
end

@defcolumn diff_confirmed(date, confirmed) _column_diff(confirmed, date, 7)
@defcolumn diff_infected(date, infected) _column_diff(infected, date, 7)
#= begin
    n = length(date)
    result = Vector{OptFloat}(missing, n)
    for i ∈ 5:(n - 3)
        if (
            date[i - 4] + Day(7) == date[i] + Day(3) == date[i + 3] && # consistency check
            !ismissing(infected[i - 4]) &&
            isfinite(infected[i - 4]) &&
            infected[i - 4] > 0
        )
            # Geometric mean of infected[i+1]/infected[i] over 7 days
            diff_rate = (infected[i + 3] / infected[i - 4])^(1 / 7)
            if !ismissing(diff_rate) && isfinite(diff_rate)
                result[i] = diff_rate * infected[i]
            end
        end
    end
    result
end
=#

@defcolumn diff_recovered(date, recovered) _column_diff(recovered, date, 7)
@defcolumn diff_deaths(date, deaths) _column_diff(deaths, date, 7)

@defcolumn diff_closed(diff_recovered, diff_deaths) diff_recovered .+ diff_deaths
@defcolumn diff_active(diff_infected, diff_closed) diff_infected .- diff_closed

@defcolumn diff2_infected(date, diff_infected) _column_diff(diff_infected, date, 5)
@defcolumn diff2_recovered(date, diff_recovered) _column_diff(diff_recovered, date, 5)
@defcolumn diff2_deaths(date, diff_deaths) _column_diff(diff_deaths, date, 5)

@defcolumn diff2_closed(diff2_recovered, diff2_deaths) diff2_recovered .+ diff2_deaths
@defcolumn diff2_active(diff2_infected, diff2_closed) diff2_infected .- diff2_closed

@defcolumn diff3_infected(date, diff2_infected) _column_diff(diff2_infected, date, 5)
#=
@defcolumn diff3_recovered(date, diff2_recovered) _column_diff(
    diff2_recovered,
    date,
    3,
)
@defcolumn diff3_deaths(date, diff2_deaths) _column_diff(diff2_deaths, date, 3)

@defcolumn diff3_closed(diff3_recovered, diff3_deaths) diff3_recovered .+ diff3_deaths
@defcolumn diff3_active(diff3_infected, diff3_closed) diff3_infected .- diff3_closed
=#
@defcolumn diff4_infected(date, diff3_infected) _column_diff(diff3_infected, date, 5)

@defcolumn diff5_infected(date, diff4_infected) _column_diff(diff4_infected, date, 5)

@defcolumn diff6_infected(date, diff5_infected) _column_diff(diff5_infected, date, 5)

const COVID_19_μ = 0.034

@defcolumn μ_closed(deaths, closed) begin
    replace(1.0 .* deaths ./ closed, NaN => missing)
end

@defcolumn infected_per_1mi(infected, estimated_population) begin
    1.0e6 .* infected ./ estimated_population
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
@defcolumn tests_per_infected(total_tests, infected) 1.0 .* total_tests ./ infected
@defcolumn tests_per_1mi(total_tests, estimated_population) begin
    1.0e6 .* total_tests ./ estimated_population
end

function copy_pop_tests_tables(data::AbstractDict)
    pop = get(data, :world_population, nothing)
    gdp_pc = get(data, :world_gdp_per_capita, nothing)
    tests = get(data, :tests_no_duplicate, nothing)
    csse = get(data, :csse, nothing)
    brasil_io = get(data, :brasil_io, nothing)
    updated = false
    if brasil_io != nothing
        if pop != nothing && !haskey(brasil_io, :world_population)
            brasil_io[:world_population] = pop
            updated = true
        end
    end
    if csse != nothing
        if pop != nothing && !haskey(csse, :world_population)
            csse[:world_population] = pop
            updated = true
        end
        if gdp_pc != nothing && !haskey(csse, :world_gdp_per_capita)
            csse[:world_gdp_per_capita] = gdp_pc
            updated = true
        end
        if tests != nothing && !haskey(csse, :tests_no_duplicate)
            csse[:tests_no_duplicate] = tests
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
        !haskey(brasil_io, :brazil) &&
        haskey(csse, :per_country_group)
    )
        brasil_io[:brazil] = csse[:per_country_group][:brazil]
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

_sum(itr) = isempty(itr) ? missing : sum(itr)
_maximum(itr) = isempty(itr) ? missing : maximum(itr)

function covid19_database(; kwargs...)
    sources = [:csse, :brasil_io, :world_population, :world_gdp_per_capita, :total_tests]
    funcs = [
        # Tables
        group_world2,
        table_tests_no_duplicate,
        copy_pop_tests_tables,
        table_world3,
        group_per_country_aux,
        group_per_country,
        group_per_region_aux,
        group_per_region,
        group_per_income_group_aux,
        group_per_income_group,
        copy_brazil_tables,
        table_brazil2,
        table_states2,
        table_cities2,
        group_per_state,
        group_per_city,
        # Columns
        column_infected,
        column_closed,
        column_active,
        column_diff_confirmed,
        column_diff_infected,
        column_diff_recovered,
        column_diff_deaths,
        column_diff_closed,
        column_diff_active,
        column_diff2_infected,
        column_diff2_recovered,
        column_diff2_deaths,
        column_diff2_closed,
        column_diff2_active,
        column_diff3_infected,
        column_diff4_infected,
        column_diff5_infected,
        # These tables needs diff columns before grouping
        group_per_date,
        table_per_date2,
        # More columns
        column_μ_closed,
        column_infected_per_1mi,
        column_deaths_per_1mi,
        column_recovered_per_1mi,
        column_closed_per_1mi,
        column_active_per_1mi,
        column_tests_per_confirmed,
        column_tests_per_1mi,
    ]
    Database{DataDict}(sources, kwargs, funcs; name = "COVID-19")
end

function covid19(;
    database = nothing,
    compute_model::Bool = true,
    pretty::Bool = true,
    do_export::Bool = true,
    export_trials = 3,
    kwargs...,
)
    db = database
    @info "Creating COVID-19 database."
    isnothing(database) && (db = covid19_database(; kwargs...))

    if hasproperty(db, :csse)
        sources = DataFrame(
            :source => [
                "CSSE at JHU (John Hopkins University)",
                "Brasil.IO",
                "World Bank Open Data",
                "Our World In Data",
            ],
            :url => [
                "https://systems.jhu.edu/",
                "https://brasil.io/home/",
                "https://data.worldbank.org/",
                "https://ourworldindata.org/",
            ],
            :data_kind => [
                "COVID-19 cases for several countries",
                "COVID-19 cases per state and city in Brazil",
                "Population and GDP per capita for several countries",
                "Number of tests of COVID-19 for several countries",
            ],
        )
        for row ∈ eachrow(db.brasil_io.estimates)
            data_kind = row.state == "total" ? "Estimate for Brazil on $(row.date)" :
                "Estimate for $(row.state), Brazil on $(row.date)"
            push!(sources, (; source = row.source, url = row.url, data_kind = data_kind))
        end

        newdata = DataDict(
            :world => DataDict(
                :per_country => db.csse.per_country_group,
                :per_region => db.csse.per_region_group,
                :per_income_group => db.csse.per_income_group_group,
                :per_date => db.csse.per_date_group,
                :total => db.csse.per_date2,
            ),
            :brazil => DataDict(
                :per_state => db.brasil_io.per_state_group,
                :per_city => db.brasil_io.per_city_group,
                :total => db.brasil_io.brazil2,
            ),
            :sources => sources,
            :auxiliar => DataDict(
                :total_tests => db.total_tests,
                :world_population => db.world_population,
                :world_gdp_per_capita => db.world_gdp_per_capita,
            ),
        )
        datadict!(db, newdata)
        br = db.brazil.total
        db.world.per_country.datadict[:brazil] = br
        @info "COVID-19 database created." summary(db)

    else
        @info "COVID-19 database imported from cache."
    end

    if compute_model || optimize
        @info "Computing SEIR model..."
        SEIRModel!(db; kwargs...)
        @info "SEIR model for COVID-19 database computed."
    end

    if do_export
        for trial ∈ 1:export_trials
            try
                @info "Exporting..."
                paths = export_data(db; pretty = pretty, kwargs...)
                break
                @info "COVID-19 database exported." paths
            catch exception
                @warn "Failed exporting data." exception trial = "$trial of $export_trials"
                sleep(1)
            end
        end
    end
    db
end
