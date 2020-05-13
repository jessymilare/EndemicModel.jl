# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

# CSSE COVID-19 Repository

const CSSE_COLUMNS = (;
    Symbol("Country/Region") => :country,
    Symbol("Province/State") => :state,
    Symbol("Lat") => :latitude,
    Symbol("Long") => :longitude,
)

function _csse_url(tablename::String)
    "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/" *
    "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_" *
    tablename *
    "_global.csv"
end

const CSSE_TABLES = (
    confirmed = _csse_url("confirmed"),
    deaths = _csse_url("deaths"),
    recovered = _csse_url("recovered"),
)

"""
Novel Coronavirus (COVID-19) Cases, provided by JHU CSSE
https://github.com/CSSEGISandData/COVID-19
https://systems.jhu.edu/research/public-health/ncov/
"""
struct CsseSource <: AbstractDataSource
    function CsseSource()
        new()
    end
end
CsseSource(::Nothing; kwargs...) = CsseSource()

DATA_SOURCES[:csse] = CsseSource

function _import_data(source::CsseSource; kwargs...)
    result = DataDict()
    for (dfsym, url) ∈ pairs(CSSE_TABLES)
        file = import_data(DownloadSource(url))
        df = csv_read(file; copycols = true)

        newcolnames = [
            if all(c -> c ∈ "0123456789/", string(colname))
                date_value = Date(string(colname), "m/d/y") + Year(2000)
                Symbol(string(date_value))
            else
                CSSE_COLUMNS[Symbol(colname)]
            end for colname ∈ names(df)
        ]
        rename!(df, newcolnames)
        df[!, :country] = correct_countries!(df.country)
        result[dfsym] = df
    end
    result
end

function import_data(source::CsseSource; kwargs...)
    @debug "Importing CSSE data."
    raw_data = _import_data(source)
    key_columns = [:country, :state, :latitude, :longitude]
    dataframes = []

    for (main_col_name, df) ∈ raw_data
        dates = []
        date_syms = []
        for colsym ∈ names(df)
            colstr = string(colsym)
            if all(c -> c ∈ "0123456789-", colstr)
                push!(dates, Date(colstr))
                push!(date_syms, colsym)
            end
        end
        datadf = combine(groupby(df, key_columns)) do subdf
            col_values = [sum(subdf[!, date_sym]) for date_sym ∈ date_syms]
            DataFrame(; date = dates, main_col_name => col_values)
        end
        push!(dataframes, datadf)
    end
    parsedata = innerjoin(dataframes...; on = key_columns ∪ [:date]) |> DataFrame
    sort!(parsedata, [:country, :state, :date])
    DataDict(:world => parsedata)
end

"""
Number of tests per country from Our World in Data.
https://ourworldindata.org/covid-testing
"""
struct TotalTestsOWID <: AbstractDataSource
    function TotalTestsOWID()
        new()
    end
end
TotalTestsOWID(::Nothing; kwargs...) = TotalTestsOWID()

DATA_SOURCES[:total_tests] = TotalTestsOWID

const TESTING_OWID_URL = "https://github.com/owid/covid-19-data/raw/master/public/data/testing/covid-testing-all-observations.csv"

function import_data(source::TotalTestsOWID; url = TESTING_OWID_URL, kwargs...)
    @debug "Importing testing data from Our World In Data."
    file = import_data(DownloadSource(url))
    col_select = [:Entity, :Date, Symbol("Cumulative total")]
    data::DataFrame = csv_read(file; copycols = true, select = col_select)
    rename!(data, [:country, :date, :total_tests])
    country_split = map(name -> split(name, " - "), data.country)
    data[!, :country] = map(row -> row[1], country_split) |> correct_countries!
    test_kinds = map(row -> row[2], country_split) |> correct_test_kind!
    insertcols!(data, 3, :test_kind => test_kinds)

    rows = collect(eachrow(data))
    nextrows = rows[2:end]
    lastrow = rows[end]
    nlastrow = DataFrame(;
        country = [lastrow.country],
        date = [today()],
        test_kind = [lastrow.test_kind],
        total_tests = [lastrow.total_tests],
    )
    push!(nextrows, nlastrow[1, :])
    growth_rate, next_growth_rate = 1.0, 1.0
    for (row, nextrow) ∈ zip(rows, nextrows)
        if row.country == nextrow.country
            dtend = nextrow.date - Day(1)
            ndays = Dates.days(nextrow.date - row.date)
            growth_rate = (nextrow.total_tests / row.total_tests)^(1.0 / ndays)
            next_growth_rate = growth_rate
        else
            dtend = today() - Day(1)
            next_growth_rate = 1.0
        end

        for (i, dt) ∈ enumerate((row.date + Day(1)):Day(1):dtend)
            tests_interpolated = round(Int, row.total_tests * growth_rate^i)
            push!(
                data,
                (;
                    country = row.country,
                    date = dt,
                    test_kind = row.test_kind,
                    total_tests = tests_interpolated,
                ),
            )
        end

        growth_rate = next_growth_rate
    end
    sort!(data, [:country, :test_kind, :date])
    data
end

# Brasil.io COVID-19 database

const BRASIL_IO_URL = "https://brasil.io/dataset/covid19/caso?format=csv"

"""
Novel Coronavirus (COVID-19) Cases in Brasil, provided by Brasil.io
https://brasil.io/dataset/covid19/caso_full/
"""
struct BrasilIo <: AbstractDataSource
    function BrasilIo()
        new()
    end
end
BrasilIo(::Nothing; kwargs...) = BrasilIo()

DATA_SOURCES[:brasil_io] = BrasilIo

const BRAZIL_POP_FACTOR = 1.01

function import_data(source::BrasilIo; url = BRASIL_IO_URL, kwargs...)
    @debug "Importing data from Brasil.io."
    col_drop = [:confirmed_per_100k_inhabitants, :death_rate, :is_last]
    file = import_data(DownloadSource(url))
    raw_data::DataFrame = csv_read(
        file;
        copycols = true,
        truestrings = ["True"],
        falsestrings = ["False"],
        drop = col_drop,
    )
    raw_data[(raw_data.deaths .=== missing), :deaths] .= 0
    sort!(raw_data, [:state, :city, :date])

    rename!(raw_data, "estimated_population_2019" => :estimated_population)
    raw_data.estimated_population =
        round.(Union{Int, Missing}, raw_data.estimated_population .* BRAZIL_POP_FACTOR)
    states = @where(raw_data, :place_type .== "state")
    cities = @where(raw_data, :place_type .== "city")

    select!(states, Not([:city, :place_type]))
    select!(cities, Not([:place_type]))

    DataDict(:states => states, :cities => cities)
end
