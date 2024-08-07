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
using LinearAlgebra
using Statistics
using Setfield
using Unitful
using Random
using Printf

import Base: merge, +, *
import Meshes: isisotropic

include("utils.jl")
include("empirical.jl")
include("variogram.jl")
include("covariance.jl")
include("nesting.jl")
include("fitting.jl")
include("sampling.jl")
include("varioplot.jl")

# temporary fix for ⋅ with missing values
# https://github.com/JuliaLang/julia/issues/40743
import LinearAlgebra: ⋅
⋅(::Missing, ::Missing) = missing

export
  # empirical variograms
  EmpiricalVariogram,
  EmpiricalVarioplane,
  DirectionalVariogram,
  PlanarVariogram,
  distance,
  estimator,

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
  sill,
  nugget,
  metricball,
  isstationary,
  isisotropic,
  structures,
  variosample,

  # theoretical covariance
  CircularCovariance,
  CubicCovariance,
  ExponentialCovariance,
  GaussianCovariance,
  MaternCovariance,
  PentasphericalCovariance,
  SineHoleCovariance,
  SphericalCovariance,

  # fitting algorithms
  VariogramFitAlgo,
  WeightedLeastSquares,

  # plotting
  varioplot,
  varioplot!

end
