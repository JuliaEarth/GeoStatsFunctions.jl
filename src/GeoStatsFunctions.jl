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

import Base: merge, +, *
import LinearAlgebra: issymmetric

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
  EmpiricalGeoStatsFunction,
  EmpiricalVariogram,
  EmpiricalTransiogram,

  # convenience functions
  DirectionalVariogram,
  DirectionalTransiogram,
  PlanarVariogram,
  PlanarTransiogram,

  # empirical surfaces
  EmpiricalGeoStatsSurface,
  EmpiricalVariogramSurface,
  EmpiricalTransiogramSurface,

  # theoretical functions
  GeoStatsFunction,
  isstationary,
  isisotropic,
  issymmetric,
  isbanded,
  metricball,
  nvariates,
  structures,
  meanlengths,
  proportions,
  sill,
  nugget,

  # theoretical variograms
  Variogram,
  GaussianVariogram,
  SphericalVariogram,
  ExponentialVariogram,
  MaternVariogram,
  CubicVariogram,
  PentaSphericalVariogram,
  SineHoleVariogram,
  CircularVariogram,
  PowerVariogram,
  NuggetEffect,

  # theoretical covariances
  Covariance,
  GaussianCovariance,
  SphericalCovariance,
  ExponentialCovariance,
  MaternCovariance,
  CubicCovariance,
  PentaSphericalCovariance,
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
  CarleTransiogram,

  # fitting algorithms
  WeightedLeastSquares,

  # plotting functions
  funplot,
  funplot!,
  surfplot,
  surfplot!

end
