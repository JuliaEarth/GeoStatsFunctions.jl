using GeoStatsFunctions
using Meshes
using GeoTables
using Unitful
using Distances
using LinearAlgebra
using CSV, DelimitedFiles
using Test, StableRNGs
import CoDa: Composition

# environment settings
datadir = joinpath(@__DIR__, "data")

# list of tests
testfiles = ["empirical.jl", "variogram.jl", "covariance.jl", "nesting.jl", "fitting.jl", "sampling.jl"]

@testset "GeoStatsFunctions.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
