# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsFunctionsMakieExt

using GeoStatsFunctions
using Meshes
using Unitful
using Distances
using TableTransforms
using LinearAlgebra
using Statistics

using GeoStatsFunctions: isinvalid

import Makie

import GeoStatsFunctions: funplot, funplot!
import GeoStatsFunctions: surfplot, surfplot!
import GeoStatsFunctions: hscatter

# source code
include("utils.jl")
include("funplot.jl")
include("surfplot.jl")
include("hscatter.jl")

# precompile workloads
include("precompile.jl")

end
