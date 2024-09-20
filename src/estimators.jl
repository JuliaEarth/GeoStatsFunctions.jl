# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    VariogramEstimator

A (robust) estimator of [`EmpiricalVariogram`](@ref).
"""
abstract type VariogramEstimator end

returntype(estim::VariogramEstimator, z₁, z₂) = typeof(formula(estim, z₁[1], z₁[2], z₂[1], z₂[2]))

include("estimators/matheron.jl")
include("estimators/cressie.jl")
