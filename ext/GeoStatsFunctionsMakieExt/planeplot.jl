# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

planeplot(f; kwargs...) = _planeplot(f; kwargs...)

include("planeplot/variogram.jl")
include("planeplot/transiogram.jl")
