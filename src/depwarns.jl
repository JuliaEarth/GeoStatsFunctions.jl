# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function EmpiricalVariogram(data::AbstractGeoTable, vars=1:(ncol(data) - 1); kwargs...)
  Base.depwarn(
    """
    `EmpiricalVariogram(data, [vars]; [options])` is deprecated.

    Use `variogram(data, [vars]; [options])` instead.
    """,
    :EmpiricalVariogram,
    force=true
  )
  variogram(data, vars; kwargs...)
end

export DirectionalVariogram

function DirectionalVariogram(dir, data::AbstractGeoTable, vars=1:(ncol(data) - 1); dtol=1e-6u"m", kwargs...)
  Base.depwarn(
    """
    `DirectionalVariogram(dir, data, [vars]; dtol=dtol, [options])` is deprecated.

    Use `variogram(data, [vars]; dir=dir, dirtol=dtol, [options])` instead.
    """,
    :DirectionalVariogram,
    force=true
  )
  variogram(data, vars; dir=dir, dirtol=dtol, kwargs...)
end

export DirectionalTransiogram

function EmpiricalTransiogram(data::AbstractGeoTable, var=1; kwargs...)
  Base.depwarn(
    """
    `EmpiricalTransiogram(data, [var]; [options])` is deprecated.

    Use `transiogram(data, [var]; [options])` instead.
    """,
    :EmpiricalTransiogram,
    force=true
  )
  transiogram(data, var; kwargs...)
end

function DirectionalTransiogram(dir, data::AbstractGeoTable, var=1; dtol=1e-6u"m", kwargs...)
  Base.depwarn(
    """
    `DirectionalTransiogram(dir, data, [var]; dtol=dtol, [options])` is deprecated.

    Use `transiogram(data, [var]; dir=dir, dirtol=dtol, [options])` instead.
    """,
    :DirectionalTransiogram,
    force=true
  )
  transiogram(data, var; dir=dir, dirtol=dtol, kwargs...)
end

function EmpiricalVariogramSurface(
  data::AbstractGeoTable,
  vars=1:(ncol(data) - 1);
  ptol=0.5u"m",
  dtol=0.5u"m",
  kwargs...
)
  Base.depwarn(
    """
    `EmpiricalVariogramSurface(data, [vars]; ptol=0.5u"m", dtol=0.5u"m", [options])` is deprecated.

    Use `variogramsurface(data, [vars]; planetol=ptol, dirtol=dtol, [options])` instead.
    """,
    :EmpiricalVariogramSurface,
    force=true
  )
  variogramsurface(data, vars; planetol=ptol, dirtol=dtol, kwargs...)
end

function EmpiricalTransiogramSurface(data::AbstractGeoTable, var=1; ptol=0.5u"m", dtol=0.5u"m", kwargs...)
  Base.depwarn(
    """
    `EmpiricalTransiogramSurface(data, [var]; ptol=0.5u"m", dtol=0.5u"m", [options])` is deprecated.

    Use `transiogramsurface(data, [var]; planetol=ptol, dirtol=dtol, [options])` instead.
    """,
    :EmpiricalTransiogramSurface,
    force=true
  )
  transiogramsurface(data, var; planetol=ptol, dirtol=dtol, kwargs...)
end
