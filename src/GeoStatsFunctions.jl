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
using NearestNeighbors: MinkowskiMetric
using OhMyThreads: tmapreduce
using DataScienceTraits
using TableTransforms
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

# empirical functions
include("empirical.jl")

# theoretical functions
include("theoretical.jl")

# misc operations
include("fitting.jl")
include("plotting.jl")

include("precompile.jl")

export
  # empirical functions
  EmpiricalVariogram,
  EmpiricalTransiogram,

  # convenience functions
  DirectionalVariogram,
  DirectionalTransiogram,
  PlanarVariogram,
  PlanarTransiogram,
  EmpiricalVarioplane,
  EmpiricalTransioplane,

  # theoretical functions
  GeoStatsFunction,
  isstationary,
  isisotropic,
  metricball,
  meanlengths,

  # theoretical variograms
  Variogram,
  GaussianVariogram,
  SphericalVariogram,
  ExponentialVariogram,
  MaternVariogram,
  CubicVariogram,
  PentasphericalVariogram,
  SineHoleVariogram,
  CircularVariogram,
  PowerVariogram,
  NuggetEffect,
  NestedVariogram,
  sill,
  nugget,
  structures,

  # theoretical covariances
  Covariance,
  GaussianCovariance,
  SphericalCovariance,
  ExponentialCovariance,
  MaternCovariance,
  CubicCovariance,
  PentasphericalCovariance,
  SineHoleCovariance,
  CircularCovariance,

  # theoretical transiograms
  Transiogram,
  LinearTransiogram,
  GaussianTransiogram,
  SphericalTransiogram,
  ExponentialTransiogram,
  MatrixExponentialTransiogram,
  PiecewiseLinearTransiogram,

  # fitting algorithms
  WeightedLeastSquares,

  # plotting functions
  varioplot,
  varioplot!,
  transioplot,
  transioplot!,
  planeplot,
  planeplot!

end
