# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# --------------------------
# CORE TYPES AND ALGORITHMS
# --------------------------

include("empirical/matrices.jl")
include("empirical/lagsearch.jl")
include("empirical/estimators.jl")

# -----------------------------
# END-USER TYPES AND FUNCTIONS
# -----------------------------

"""
    EmpiricalGeoStatsFunction

An empirical geostatistical function.
"""
abstract type EmpiricalGeoStatsFunction end

"""
    issymmetric(f)

Tells whether or not the empirical geostatistical function `f` is symmetric.
"""
issymmetric(f::EmpiricalGeoStatsFunction) = issymmetric(typeof(f))

"""
    nvariables(f)

Return the number of (co)variables of the empirical geostatistical function `f`.
"""
nvariables(f::EmpiricalGeoStatsFunction) = nvariables(typeof(f))

"""
    variables(f)

Return the names of the (co)variables of the empirical geostatistical function `f`.
"""
variables(f::EmpiricalGeoStatsFunction) = defaultvariables(nvariables(f))

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, f::EmpiricalGeoStatsFunction) = print(io, nameof(typeof(f)))

function Base.show(io::IO, f::EmpiricalGeoStatsFunction)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, f)
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
  print(ioctx, "├─ variables: ")
  _printlnvec(ioctx, f.variables, 3)
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

Tells whether or not the empirical geostatistical surface `f` is symmetric.
"""
issymmetric(f::EmpiricalGeoStatsSurface) = issymmetric(typeof(f))

"""
    nvariables(f)

Return the number of (co)variables of the empirical geostatistical surface `f`.
"""
nvariables(f::EmpiricalGeoStatsSurface) = nvariables(typeof(f))

"""
    variables(f)

Return the names of the (co)variables of the empirical geostatistical surface `f`.
"""
variables(f::EmpiricalGeoStatsSurface) = defaultvariables(nvariables(f))

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, f::EmpiricalGeoStatsSurface) = print(io, nameof(typeof(f)))

function Base.show(io::IO, f::EmpiricalGeoStatsSurface)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, f)
end

function Base.show(io::IO, ::MIME"text/plain", f::EmpiricalGeoStatsSurface)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, f)
  println(ioctx)
  println(ioctx, "├─ nangs: ", length(f.θs))
  print(ioctx, "└─ nlags: ", length(f.rs))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("empirical/variosurf.jl")
include("empirical/transiosurf.jl")
