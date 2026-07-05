# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    fit(F, f, w=inv; parameters..., maxparameters...)

Fit theoretical geostatistical function of type `F` to empirical function `f`.
The weighting function `w` is optional and maps lag distances to importance weights.

Theoretical `parameters` of `F` like `range`, `sill` and `nugget` as well as their
maximum values like `maxrange`, `maxsill` and `maxnugget` can be fixed as keyword
arguments. They are converted into optimization constraints in the weighted least
squares fitting algorithm.

The type `F` can be abstract like `Variogram` or `Transiogram`, in which case
all fittable (e.g., stationary) subtypes are fitted and the one with minimum
error is returned. See [`fitany`](@ref) for custom lists of types.

    fit([F₁, F₂, ..., Fₙ], f, w=inv; parameters..., maxparameters...)

Alternatively, fit linear model of coregionalization (LMC) with structures of type
`F₁, F₂, ..., Fₙ`. The result is a [`CompositeFunction`](@ref) that can be written as
`f = Bₒfₒ + B₁f₁ + B₂f₂ + ⋯ + Bₙfₙ` where `fₒ` is a [`NuggetEffect`](@ref) and
`Bₒ, B₁, B₂, ..., Bₙ` are positive semidefinite coefficient matrices.

## Examples

```julia
julia> fit(SphericalVariogram, g)
julia> fit(ExponentialVariogram, g, range=1.0)
julia> fit(GaussianVariogram, g, sill=1.0, nugget=0.1)
julia> fit(ExponentialVariogram, g, maxsill=1.0)
julia> fit(SphericalVariogram, g, h -> 1/h^2)
julia> fit(Variogram, g, range=1.0)
julia> fit(Transiogram, t, h -> exp(-h))
julia> fit([SphericalVariogram, ExponentialVariogram], g)
```
"""
fit(F::Type{<:GeoStatsFunction}, f::EmpiricalGeoStatsFunction, w=inv; kwargs...) = _fit(F, f, w; kwargs...) |> first

fit(Fs::AbstractVector, f::EmpiricalGeoStatsFunction, w=inv; kwargs...) = _fitlmc(Fs, f, w; kwargs...) |> first

function fit(::Type{Variogram}, g::EmpiricalVariogram, w=inv; kwargs...)
  Gs = (
    CircularVariogram,
    CubicVariogram,
    ExponentialVariogram,
    GaussianVariogram,
    MaternVariogram,
    PentaSphericalVariogram,
    SineHoleVariogram,
    SphericalVariogram
  )
  fitany(Gs, g, w; kwargs...)
end

function fit(::Type{Transiogram}, t::EmpiricalTransiogram, w=inv; kwargs...)
  Ts =
    (ExponentialTransiogram, GaussianTransiogram, LinearTransiogram, MatrixExponentialTransiogram, SphericalTransiogram)
  fitany(Ts, t, w; kwargs...)
end

"""
    fitany([F₁, F₂, ..., Fₙ], f, w=inv; parameters..., maxparameters...)

Fit theoretical geostatistical functions of types `F₁, F₂, ..., Fₙ`
to empirical function `f` and return the one with minimum error.

Please check [`fit`](@ref) for more detailed documentation on the
arguments and keyword arguments.

## Examples

```julia
julia> fitany([SphericalVariogram, ExponentialVariogram], g)
```
"""
function fitany(Fs, f::EmpiricalGeoStatsFunction, w=inv; kwargs...)
  sols = [_fit(F, f, w; kwargs...) for F in Fs]
  f, _ = argmin(last, sols)
  f
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("fitting/variograms.jl")
include("fitting/transiograms.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

_ustrip(u, x) = x
_ustrip(u, x::Quantity) = ustrip(u, x)
_ustrip(u, X::AbstractMatrix{<:Quantity}) = ustrip.(u, X)

_weights(w, x) = map(xᵢ -> ustrip(w(xᵢ)), x)

function _coefmat(θ, k)
  L = _lowertrimat(θ, k)
  Symmetric(L * transpose(L))
end

function _lowertrimat(θ, k)
  T = eltype(θ)
  p(i, j) = i + (j - 1) * k - (j - 1) * j ÷ 2
  A = SMatrix{k,k}(i ≥ j ? θ[p(i, j)] : zero(T) for i in 1:k, j in 1:k)
  LowerTriangular(A)
end

function _coefvec(A)
  L = cholesky(A).L
  k = size(L, 1)
  [L[i, j] for j in 1:k for i in j:k]
end

function _posmat(A)
  λ, V = eigen((A + transpose(A)) / 2)
  δ = eltype(λ)(1e-8)
  P = V * Diagonal(max.(λ, δ)) * transpose(V)
  Symmetric((P + transpose(P)) / 2)
end

function _optimize(J, θₗ, θᵤ, θₒ)
  s = Optim.optimize(J, θₗ, θᵤ, θₒ, LBFGSB())
  ϵ = Optim.minimum(s)
  θ = Optim.minimizer(s)
  θ, ϵ
end
