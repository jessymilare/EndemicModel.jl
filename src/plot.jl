# This file is part of EndemicModel project, which is licensed under BDS-3-Clause.
# See file LICENSE.md for more information.

const FACTOR_LABEL_MAP = OrderedDict(
    1e8 => N_"hundreds of millions",
    1e7 => N_"tens of millions",
    1e6 => N_"millions",
    1e5 => N_"hundreds of thousands",
    1e4 => N_"tens of thousands",
    1e3 => N_"thousands",
    1e2 => N_"hundreds",
    1e1 => N_"tens",
    1.0 => "",
)

function _get_yscale(max_yvalue)
    for (yscale, ylabel) ∈ FACTOR_LABEL_MAP
        if max_yvalue / yscale >= 2.0
            return (yscale, gettext(ylabel))
        end
    end
    (1.0, "")
end

function _get_plot_title(df)
    city = "city" ∈ names(df) ? df[1, :city] : missing
    state = "state" ∈ names(df) ? df[1, :state] : missing
    country = "country" ∈ names(df) ? df[1, :country] : missing
    location = skipmissing([city, state, country])
    join(location, ", ")
end

plot(args...; kwargs...) = Plots.plot(args...; kwargs...)

function plot(
    df::DataFrame;
    columns = option(:plot_columns),
    plot_begin::Union{Integer,Date,Period} = 1,
    plot_end::Union{Integer,Date,Period} = nrow(df),
    labels = nothing,
    title = _get_plot_title(df),
    date_format = option(:plot_date_format),
    new_window::Bool = true,
    ylabel = _"People",
    yscale = nothing,
    legend = length(columns) == 1 ? false : :topleft,
    left_margin = 1mm,
    kwargs...,
)
    maybe_load_language()

    if plot_begin isa Date
        ini = findfirst(isequal(plot_begin), df.date)
        if isnothing(ini)
            throw(ArgumentError("plot_begin not found in DataFrame: $(plot_begin)"))
        end
    elseif plot_begin isa Period
        ini = Dates.days(plot_begin) + 1
    elseif plot_begin <= 0
        ini = max(1, nrow(df) + plot_begin)
    else
        ini = min(plot_begin, nrow(df))
    end

    if plot_end isa Date
        iend = findfirst(isequal(plot_end), df.date)
        if isnothing(ini)
            throw(ArgumentError("plot_end not found in DataFrame: $(plot_end)"))
        end
    elseif plot_end isa Period
        iend = nrow(df) - Dates.days(plot_end)
    elseif plot_end <= 0
        iend = max(1, nrow(df) + plot_end)
    else
        iend = min(plot_end, nrow(df))
    end
    @argcheck ini <= iend

    if isa(columns, String) || isa(columns, Symbol)
        columns = [columns]
    end

    df = df[ini:iend, :]
    columns = intersect(string.(columns), names(df))
    if isnothing(labels)
        labels = gettext.(string.(prettify.(Symbol.(columns))))
    end
    max_yvalue = 0
    for cname ∈ columns
        max_yvalue = max(maximum(skipmissing(df[!, cname])), max_yvalue)
    end
    numpeople = if hasproperty(df, :estimated_population)
        df.estimated_population[1]
    else
        2 * max_yvalue
    end
    if isnothing(yscale)
        (yscale, scalestr) = _get_yscale(max_yvalue)
        !isempty(scalestr) && (ylabel *= " ($(scalestr))")
    end

    X = Dates.format.(df.date, date_format)

    left_margin isa Real && (left_margin *= mm)

    win = Plots.plot(
        X,
        df[!, columns[1]] ./ yscale,
        xrotation = 45,
        xticks = 15,
        label = labels[1],
        ylabel = ylabel,
        title = title,
        legend = legend,
        left_margin = left_margin,
    )
    for (yn, label) ∈ zip(columns[2:end], labels[2:end])
        Plots.plot!(win, X, df[!, yn] ./ yscale; label = label)
    end
    new_window && gui(win)
    win
end

function model_combined_data(
    model::AbstractEndemicModel;
    columns = option(:plot_columns),
    plot_period = Day(14),
    kwargs...,
)
    columns = string.(columns)

    realdf = realdata(model)
    model_initial_date = default_kwarg(model, :initial_date, today() - Day(15))
    plot_initial_date = realdf.date[end] - plot_period
    ndays = Dates.days(model_initial_date - plot_initial_date)

    @debug "Combining data" model_initial_date plot_initial_date ndays

    realdf = realdf[(end-Dates.days(plot_period)):end, :]
    if ndays > 1
        model = model_step(model, -ndays; initial_date = plot_initial_date)
    end
    modeldf = @where(
        modeldata(model),
        (:date .<= realdf.date[end]) .& (:date .>= plot_initial_date)
    )

    columns = intersect(columns, names(modeldf), names(realdf))
    rcolumns = [Symbol(sym, " (real)") for sym ∈ columns]
    mcolumns = [Symbol(sym, " (model)") for sym ∈ columns]

    realdf = select(realdf, :date, (columns .=> rcolumns)...)
    modeldf = select(modeldf, :date, (columns .=> mcolumns)...)

    @debug "Model dataframes computed" _debuginfo(realdf) _debuginfo(modeldf)
    (leftjoin(realdf, modeldf; on = :date), [rcolumns; mcolumns])
end

function plot(
    model::AbstractEndemicModel,
    kind = :combined_data;
    columns = option(:plot_columns),
    new_window::Bool = true,
    title = nothing,
    kwargs...,
)
    kind = Symbol(kind)
    if kind == :combined_data
        (df, columns) = model_combined_data(model; columns = columns, kwargs...)
    elseif kind == :realdata
        df = realdata(model)
    elseif kind == :modeldata
        df = modeldata(model)
    else
        throw(ArgumentError("Unrecognized plot kind: $(kind)"))
    end
    isnothing(title) && (title = _get_plot_title(realdata(model)))

    plot(df; title = title, new_window = new_window, columns = columns, kwargs...)
end
