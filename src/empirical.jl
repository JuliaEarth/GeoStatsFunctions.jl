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

"""
    distance(f)

Distance used to compute the empirical function `f`.
"""
distance(f::EmpiricalFunction) = f.distance

"""
    estimator(f)

Estimator used to compute the empirical function `f`.
"""
estimator(f::EmpiricalFunction) = f.estimator

include("empirical/variogram.jl")
include("empirical/varioplane.jl")
include("empirical/transiogram.jl")
