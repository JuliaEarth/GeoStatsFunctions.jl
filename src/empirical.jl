# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# --------------------------
# CORE TYPES AND ALGORITHMS
# --------------------------

include("empirical/matrices.jl")
include("empirical/estimators.jl")
include("empirical/algorithms.jl")
include("empirical/estimalgo.jl")

# -----------------------------
# END-USER TYPES AND FUNCTIONS
# -----------------------------

"""
    EmpiricalGeoStatsFunction

An empirical geostatistical function.
"""
abstract type EmpiricalGeoStatsFunction end

"""
    nvariates(f)

Return the number of (co)variates of the empirical geostatistical function `f`.
"""
nvariates(f::EmpiricalGeoStatsFunction) = nvariates(typeof(f))

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, f::EmpiricalGeoStatsFunction) = print(io, nameof(typeof(f)))

function Base.show(io::IO, f::EmpiricalGeoStatsFunction)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, f)
  print(ioctx, "(")
  print(ioctx, "distance: ", f.distance)
  print(ioctx, ", estimator: ", f.estimator)
  print(ioctx, ", npairs: ", sum(f.counts))
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", f::EmpiricalGeoStatsFunction)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, f)
  println(ioctx)
  print(ioctx, "├─ abscissas: ")
  _printlnvec(ioctx, f.abscissas, 3)
  print(ioctx, "├─ ordinates: ")
  _printlnvec(ioctx, f.ordinates, 3)
  println(ioctx, "├─ distance: ", f.distance)
  println(ioctx, "├─ estimator: ", f.estimator)
  print(ioctx, "└─ npairs: ", sum(f.counts))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("empirical/variogram.jl")
include("empirical/transiogram.jl")

# ---------------------
# HIGH-ORDER FUNCTIONS
# ---------------------

"""
    EmpiricalGeoStatsSurface

An empirical geostatistical surface.
"""
abstract type EmpiricalGeoStatsSurface end

"""
    issymmetric(f)

Tell whether or not the empirical geostatistical surface `f` is symmetric.
"""
issymmetric(f::EmpiricalGeoStatsSurface) = issymmetric(typeof(f))

"""
    nvariates(f)

Return the number of (co)variates of the empirical geostatistical surface `f`.
"""
nvariates(f::EmpiricalGeoStatsSurface) = nvariates(typeof(f))

include("empirical/variosurf.jl")
include("empirical/transiosurf.jl")
