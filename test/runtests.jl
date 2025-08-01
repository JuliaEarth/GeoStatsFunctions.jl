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
using CoDa: Composition

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = [
  # empirical
  "empirical/matrices.jl",
  "empirical/variogram.jl",
  "empirical/transiogram.jl",
  "empirical/variosurf.jl",
  "empirical/transiosurf.jl",

  # theoretical
  "theoretical/variogram.jl",
  "theoretical/covariance.jl",
  "theoretical/transiogram.jl",
  "theoretical/composite.jl",

  # misc operations
  "sampling.jl",
  "fitting.jl"
]

@testset "GeoStatsFunctions.jl" begin
  for testfile in testfiles
    println("Testing $testfile...")
    include(testfile)
  end
end
