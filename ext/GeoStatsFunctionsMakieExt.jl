# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsFunctionsMakieExt

using GeoStatsFunctions
using Meshes
using Unitful
using LinearAlgebra
using Statistics

import Makie

import GeoStatsFunctions: funplot, funplot!
import GeoStatsFunctions: surfplot, surfplot!

# source code
include("utils.jl")
include("funplot.jl")
include("surfplot.jl")

# precompile workloads
include("precompile.jl")

end
