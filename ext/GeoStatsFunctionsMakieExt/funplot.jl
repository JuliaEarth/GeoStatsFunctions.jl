# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

funplot(f; kwargs...) = _funplot(f; kwargs...)

include("funplot/variogram.jl")
include("funplot/covariance.jl")
include("funplot/transiogram.jl")
