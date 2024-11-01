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
    EmpiricalFunction

An empirical function estimated from data.
"""
abstract type EmpiricalFunction end

# -----------
# IO METHODS
# -----------

Base.summary(io::IO, γ::EmpiricalFunction) = print(io, nameof(typeof(γ)))

function Base.show(io::IO, γ::EmpiricalFunction)
  ioctx = IOContext(io, :compact => true)
  summary(ioctx, γ)
  print(ioctx, "(")
  print(ioctx, "distance: ", γ.distance)
  print(ioctx, ", estimator: ", γ.estimator)
  print(ioctx, ", npairs: ", sum(γ.counts))
  print(ioctx, ")")
end

function Base.show(io::IO, ::MIME"text/plain", γ::EmpiricalFunction)
  ioctx = IOContext(io, :compact => true, :limit => true)
  summary(ioctx, γ)
  println(ioctx)
  print(ioctx, "├─ abscissas: ")
  _printlnvec(ioctx, γ.abscissas, 3)
  print(ioctx, "├─ ordinates: ")
  _printlnvec(ioctx, γ.ordinates, 3)
  println(ioctx, "├─ distance: ", γ.distance)
  println(ioctx, "├─ estimator: ", γ.estimator)
  print(ioctx, "└─ npairs: ", sum(γ.counts))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("empirical/variogram.jl")
include("empirical/transiogram.jl")

# ---------------------
# HIGH-ORDER FUNCTIONS
# ---------------------

include("empirical/varioplane.jl")
include("empirical/transioplane.jl")
