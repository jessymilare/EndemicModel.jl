# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

# World Bank Population database

const WORLD_POP_WB_URL = "http://api.worldbank.org/v2/en/indicator/SP.POP.TOTL?downloadformat=csv"
const WORLD_POP_YEAR_0 = 2000

"""
Total population of each country, extracted from https://data.worldbank.org/
"""
struct WorldPopulationWB <: AbstractDataSource
    function WorldPopulationWB()
        new()
    end
end
WorldPopulationWB(::Nothing; kwargs...) = WorldPopulationWB()

DATA_SOURCES[:world_population] = WorldPopulationWB

function import_data(
    source::WorldPopulationWB;
    update_period = Year(1),
    url = WORLD_POP_WB_URL,
    kwargs...,
)
    @debug "Importing population data from World Bank."
    curyear = year(today())
    years = [string(y) for y ∈ WORLD_POP_YEAR_0:curyear]
    keycols = ["Country Name", "Country Code"]
    col_select = [keycols; years]

    file = import_data(DownloadSource(url); update_period = update_period, kwargs...)
    global r = ZipFile.Reader(file)
    zipfiles = r.files
    file = zipfiles[findfirst(f -> f.name[1:10] == "API_SP.POP", zipfiles)]

    data = csv_read(file; copycols = true, header = 5, select = col_select)
    close(file)
    global r = nothing

    rename!(data, "Country Name" => :country, "Country Code" => :country_code)
    data[!, :country] = correct_countries!(data.country)
    sort!(data, :country)

    ncols = ncol(data)
    for nr ∈ ncols:-1:1
        if !all(ismissing.(data[!, nr]))
            if nr < ncols
                select!(data, Not((nr + 1):ncols))
            end
            break
        end
    end
    years = intersect(Symbol.(years), Symbol.(names(data)))
    fyear = years[1]
    lyear = years[end]

    delete!(data, data[!, fyear] .=== missing)
    myear = years[findlast(y -> all(data[!, y] .!== missing), years)]

    fyear_int = parse(Int, string(fyear))
    lyear_int = parse(Int, string(lyear))
    myear_int = parse(Int, string(myear))

    n = (lyear_int - fyear_int)
    avg_pop_factor = (data[!, lyear] ./ data[!, fyear]) .^ (1 / n)
    navg_factor = Float64[
        factor === missing ? (data[i, myear] / data[i, fyear])^(1 / n) : factor
        for (i, factor) ∈ enumerate(avg_pop_factor)
    ]
    insertcols!(data, ncol(data) + 1, :avg_pop_factor => navg_factor)

    n = (curyear - lyear_int)
    pop_est = (data[!, lyear] .* (data.avg_pop_factor .^ n))
    pop_est = [
        ismissing(val) ? data[i, myear] * (data[i, :avg_pop_factor] .^ n) : val
        for (i, val) ∈ enumerate(pop_est)
    ]
    pop_est = pop_est .|> round .|> Int
    insertcols!(data, ncol(data) + 1, :estimated_population => pop_est)

    select(data, :country, :country_code, :avg_pop_factor, :estimated_population)
end

# World Bank GDP per capita database

const WORLD_GDP_PER_CAPITA_WB_URL = "http://api.worldbank.org/v2/en/indicator/NY.GDP.PCAP.CD?downloadformat=csv"

const WORLD_GDP_PER_CAPITA_YEAR_0 = 2000

"""
GDP per capita of each country, extracted from https://data.worldbank.org/
"""
struct WorldGDPPerCapitaWB <: AbstractDataSource
    function WorldGDPPerCapitaWB()
        new()
    end
end
WorldGDPPerCapitaWB(::Nothing; kwargs...) = WorldGDPPerCapitaWB()

DATA_SOURCES[:world_gdp_per_capita] = WorldGDPPerCapitaWB

function import_data(
    source::WorldGDPPerCapitaWB;
    update_period = Year(1),
    url = WORLD_GDP_PER_CAPITA_WB_URL,
    kwargs...,
)
    @debug "Importing GDP per capita data from World Bank."
    curyear = year(today())
    years = [string(y) for y ∈ WORLD_GDP_PER_CAPITA_YEAR_0:curyear]
    keycols = ["Country Name", "Country Code"]
    col_select = [keycols; years]

    file = import_data(DownloadSource(url); update_period = update_period, kwargs...)
    global r = ZipFile.Reader(file)
    zipfiles = r.files
    file = zipfiles[findfirst(f -> f.name[1:15] == "API_NY.GDP.PCAP", zipfiles)]

    data = csv_read(file; copycols = true, header = 5, select = col_select)
    close(file)
    global r = nothing

    rename!(data, "Country Name" => :country, "Country Code" => :country_code)
    data[!, :country] = correct_countries!(data.country)
    sort!(data, :country)

    ncols = ncol(data)
    for nr ∈ ncols:-1:1
        if !all(ismissing.(data[!, nr]))
            if nr < ncols
                select!(data, Not((nr + 1):ncols))
            end
            break
        end
    end

    years = intersect(Symbol.(years), Symbol.(names(data)))
    fyear = years[1]
    lyear = years[end]

    delete!(data, data[!, fyear] .=== missing)
    myear = years[findlast(y -> all(data[!, y] .!== missing), years)]

    fyear_int = parse(Int, string(fyear))
    lyear_int = parse(Int, string(lyear))
    myear_int = parse(Int, string(myear))

    n = (lyear_int - fyear_int)
    avg_gdp_pc_factor = (data[!, lyear] ./ data[!, fyear]) .^ (1 / n)
    navg_factor = Float64[
        factor === missing ? (data[i, myear] / data[i, fyear])^(1 / n) : factor
        for (i, factor) ∈ enumerate(avg_gdp_pc_factor)
    ]
    insertcols!(data, ncol(data) + 1, :avg_gdp_pc_factor => navg_factor)

    n = (curyear - lyear_int)
    gdp_pc_est = data[!, lyear]
    gdp_pc_est =
        [ismissing(val) ? data[i, myear] : val for (i, val) ∈ enumerate(gdp_pc_est)]
    gdp_pc_est = gdp_pc_est .|> round .|> Int
    insertcols!(data, ncol(data) + 1, :gdp_per_capita => gdp_pc_est)

    select(data, :country, :country_code, :avg_gdp_pc_factor, :gdp_per_capita)
end
