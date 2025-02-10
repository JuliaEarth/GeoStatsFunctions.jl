# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsFunctionsMakieExt

using GeoStatsFunctions
using Meshes
using Unitful
using LinearAlgebra

import Makie

import GeoStatsFunctions: funplot, planeplot

include("funplot.jl")
include("planeplot.jl")
include("utils.jl")

end
