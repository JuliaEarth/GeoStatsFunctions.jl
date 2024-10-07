# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsFunctionsMakieExt

using GeoStatsFunctions
using Unitful

import Makie

import GeoStatsFunctions: varioplot, varioplot!
import GeoStatsFunctions: transioplot, transioplot!

include("varioplot.jl")
include("transioplot.jl")
include("utils.jl")

end
