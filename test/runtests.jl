using GeoStatsFunctions
using CoordRefSystems
using Meshes
using GeoTables
using Unitful
using Distances
using StaticArrays
using LinearAlgebra
using CSV, DelimitedFiles
using Test, StableRNGs
import CoDa: Composition

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  # empirical
  "empirical/matrices.jl",
  "empirical/variogram.jl",
  "empirical/varioplane.jl",
  "empirical/transiogram.jl",
  "empirical/transioplane.jl",

  # theoretical
  "theoretical/variogram.jl",
  "theoretical/covariance.jl",
  "theoretical/transiogram.jl",
  "theoretical/sampling.jl",

  # misc operations
  "fitting.jl"
]

@testset "GeoStatsFunctions.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
