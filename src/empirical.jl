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

include("empirical/varioplane.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

function _printlnvec(io, vec, n)
  _printvec(io, vec, n)
  println(io)
end

function _printvec(io, vec::AbstractArray, n)
  print(io, "[")
  if length(vec) > 2n
    k = n - 1
    join(io, vec[begin:(begin + k)], ", ")
    print(io, ", ..., ")
    join(io, vec[(end - k):end], ", ")
  else
    join(io, vec, ", ")
  end
  print(io, "]")
end

function _printvec(io, vec::AbstractArray{<:AbstractArray}, n)
  len = length(vec)
  println(io)
  if len > 2n
    for i in 1:n
      print(io, "│  ├─ ")
      _printlnvec(io, vec[i], n)
    end
    println(io, "│  ⋮")
    for i in (len - n + 1):(len - 1)
      print(io, "│  ├─ ")
      _printlnvec(io, vec[i], n)
    end
  else
    for i in 1:(len - 1)
      print(io, "│  ├─ ")
      _printlnvec(io, vec[i], n)
    end
  end
  print(io, "│  └─ ")
  _printvec(io, vec[len], n)
end
