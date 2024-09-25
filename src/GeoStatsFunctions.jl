# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module GeoStatsFunctions

using Meshes
using GeoTables

using Optim
using Tables
using Distances
using Bessels: gamma, besselk
using InteractiveUtils: subtypes
using NearestNeighbors: MinkowskiMetric
using OhMyThreads: tmapreduce
using DataScienceTraits
using CategoricalArrays
using StaticArrays
using LinearAlgebra
using Statistics
using Setfield
using Unitful
using Random
using Printf

import Base: merge, +, *
import Meshes: isisotropic

# temporary fix for ⋅ with missing values
# https://github.com/JuliaLang/julia/issues/40743
import LinearAlgebra: ⋅
⋅(::Missing, ::Missing) = missing

# utilities
include("utils.jl")

# relevant matrices
include("matrices.jl")

# empirical models
include("empirical.jl")

# theoretical models
include("theoretical.jl")

# operations
include("fitting.jl")

# plotting
include("varioplot.jl")

include("precompile.jl")

export
  # empirical variograms
  EmpiricalVariogram,
  EmpiricalVarioplane,
  DirectionalVariogram,
  PlanarVariogram,
  distance,
  estimator,

  # theoretical function
  GeoStatsFunction,
  isisotropic,
  metricball,

  # theoretical variograms
  Variogram,
  NuggetEffect,
  GaussianVariogram,
  ExponentialVariogram,
  MaternVariogram,
  SphericalVariogram,
  CubicVariogram,
  PentasphericalVariogram,
  PowerVariogram,
  SineHoleVariogram,
  CircularVariogram,
  NestedVariogram,
  isstationary,
  sill,
  nugget,
  structures,

  # theoretical covariance
  Covariance,
  CircularCovariance,
  CubicCovariance,
  ExponentialCovariance,
  GaussianCovariance,
  MaternCovariance,
  PentasphericalCovariance,
  SineHoleCovariance,
  SphericalCovariance,

  # theoretical transiogram
  Transiogram,
  ExponentialTransiogram,
  ratematrix,

  # fitting algorithms
  VariogramFitAlgo,
  WeightedLeastSquares,

  # plotting
  varioplot,
  varioplot!

end
