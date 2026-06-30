# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function hscatter(
  data,
  vars;
  # lag distance in length units
  lag=0.0u"m",
  # tolerance for lag distance
  tol=0.1u"m",
  # distance from Distances.jl
  distance=Euclidean(),
  # size of points
  size=2,
  # color of points
  color=:black,
  # transparency of points
  alpha=1.0,
  # color of regression line
  rcolor=:salmon,
  # color of identity line
  icolor=:black,
  # color of center lines
  ccolor=:teal
)
  # selected variables
  sdata = data |> Select(vars)

  # pairs of variables
  svars = setdiff(names(sdata), ["geometry"])
  pairs = [(var₁, var₂) for var₁ in svars, var₂ in svars]

  n = length(svars)
  fig = Makie.Figure()
  for i in 1:n, j in 1:n
    i < j && continue
    # retrieve variable names
    var₁, var₂ = pairs[j, i]

    # initialize axis
    ax = Makie.Axis(fig[i, j])
    ax.aspect = Makie.AxisAspect(1)
    ax.xlabel = var₁
    ax.ylabel = var₂
    i < n && Makie.hidexdecorations!(ax, grid=false)
    j > 1 && Makie.hideydecorations!(ax, grid=false)

    # compute h-scatter coordinates
    z₁, z₂ = _hscatter(sdata, var₁, var₂, aslen(lag), aslen(tol), distance)

    # skip empty lag distances
    isnothing(z₁) && continue

    # compute regression line and identity line
    z̄₁, z̄₂ = mean(z₁), mean(z₂)
    Z₁ = [z₁ ones(length(z₁))]
    ẑ₂ = Z₁ * (Z₁ \ z₂)
    a, b = extrema([extrema(z₁)..., extrema(z₂)...])
    minmax = [(a, a), (b, b)]

    # plot h-scatter points
    Makie.scatter!(ax, z₁, z₂, color=color, alpha=alpha, markersize=size)

    # plot regression line
    Makie.lines!(ax, z₁, ẑ₂, color=rcolor)

    # plot identity line
    Makie.lines!(ax, minmax, color=icolor)

    # plot center lines
    Makie.vlines!(ax, z̄₁, color=ccolor)
    Makie.hlines!(ax, z̄₂, color=ccolor)

    # plot center point
    Makie.scatter!(ax, z̄₁, z̄₂, color=ccolor, marker=:rect, markersize=16)
  end

  for i in 1:n
    axes = [c for c in Makie.contents(fig[i, :]) if c isa Makie.Axis]
    Makie.linkyaxes!(axes...)
  end
  for j in 1:n
    axes = [c for c in Makie.contents(fig[:, j]) if c isa Makie.Axis]
    Makie.linkxaxes!(axes...)
  end

  fig
end

function _hscatter(data, var₁, var₂, lag, tol, distance)
  # lookup valid data
  ind₁ = findall(!isinvalid, data[:, var₁])
  ind₂ = findall(!isinvalid, data[:, var₂])
  val₁ = data[ind₁, :]
  val₂ = data[ind₂, :]

  # subsample to avoid long-waiting times
  nmax = 4000
  sub₁ = val₁ |> Sample(min(nrow(val₁), nmax), replace=false)
  sub₂ = val₂ |> Sample(min(nrow(val₂), nmax), replace=false)
  dom₁ = domain(sub₁)
  dom₂ = domain(sub₂)

  # extract coordinates and values
  c₁ = [centroid(dom₁, i) for i in 1:nelements(dom₁)]
  c₂ = [centroid(dom₂, i) for i in 1:nelements(dom₂)]
  z₁ = sub₁[:, var₁]
  z₂ = sub₂[:, var₂]

  # compute pairwise distance
  m, n = length(z₁), length(z₂)
  pairs = [(i, j) for j in 1:n for i in j:m]
  dists = [evaluate(distance, c₁[i], c₂[j]) for (i, j) in pairs]

  # find indices with given lag
  match = findall(abs.(dists .- lag) .< tol)

  # return sentinel values if no match found
  isempty(match) && return (nothing, nothing)

  # h-scatter coordinates
  ij = view(pairs, match)
  z₁[first.(ij)], z₂[last.(ij)]
end
