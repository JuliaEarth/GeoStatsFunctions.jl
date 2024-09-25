# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeoStatsFunction

A theoretical geostatistical function (e.g. variogram, covariance).
"""
abstract type GeoStatsFunction end

include("theoretical/variogram.jl")
include("theoretical/covariance.jl")
include("theoretical/transiogram.jl")
