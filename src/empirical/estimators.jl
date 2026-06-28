# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FunctionEstimator

An estimator of geostatistical functions.
"""
abstract type FunctionEstimator end

# ---------------------
# VARIOGRAM ESTIMATORS
# ---------------------

"""
    VariogramEstimator

An estimator of variogram functions.
"""
abstract type VariogramEstimator <: FunctionEstimator end

returntype(e::VariogramEstimator, z₁, z₂) = typeof(accumterm(e, z₁[1], z₁[2], z₂[1], z₂[2]))

include("estimators/matheron.jl")
include("estimators/cressie.jl")

# -----------------------
# TRANSIOGRAM ESTIMATORS
# -----------------------

"""
    TransiogramEstimator

An estimator of transiogram functions.
"""
abstract type TransiogramEstimator <: FunctionEstimator end

include("estimators/carle.jl")
