# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Variogram

A theoretical variogram function (e.g. Gaussian variogram).
"""
abstract type Variogram <: GeoStatsFunction end

"""
    sill(γ)

Return the sill of the variogram `γ`.
"""
sill(γ::Variogram) = γ.sill

"""
    nugget(γ)

Return the nugget of the variogram `γ`.
"""
nugget(γ::Variogram) = γ.nugget

"""
    structures(γ)

Return the individual structures of a (possibly composite)
variogram as a tuple. The structures are the total nugget,
and the coefficients (or contributions) for for the remaining
non-trivial structures after normalization (i.e. sill=1, nugget=0).

## Examples

```julia
γ₁ = GaussianVariogram(nugget=1, sill=2)
γ₂ = SphericalVariogram(nugget=2, sill=3)

structures(2γ₁ + 3γ₂)
```
"""
function structures(γ::Variogram)
  cₒ = nugget(γ)
  c = sill(γ) - nugget(γ)
  γ = @set γ.sill = oneunit(c)
  γ = @set γ.nugget = zero(c)
  cₒ, (c,), (γ,)
end

"""
    scale(γ, s)

Scale metric ball of variogram `γ` with strictly
positive scaling factor `s`.
"""
function scale(γ::Variogram, s::Real)
  V = constructor(γ)
  V(s * metricball(γ); sill=sill(γ), nugget=nugget(γ))
end

# leverage symmetry of variograms
function pairwise!(Γ, γ::Variogram, domain)
  _, (_, n, k) = matrixparams(γ, domain, domain)
  @inbounds for j in 1:n
    gⱼ = domain[j]
    sⱼ = _sample(γ, gⱼ)
    # lower triangular entries
    for i in (j + 1):n
      gᵢ = domain[i]
      sᵢ = _sample(γ, gᵢ)
      Γᵢⱼ = ustrip.(mean(γ(pᵢ, pⱼ) for pᵢ in sᵢ, pⱼ in sⱼ))
      Γ[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k)] .= Γᵢⱼ
    end
    # diagonal entries
    Γᵢⱼ = ustrip.(mean(γ(pⱼ, pⱼ) for pⱼ in sⱼ, pⱼ in sⱼ))
    Γ[((j - 1) * k + 1):(j * k), ((j - 1) * k + 1):(j * k)] .= Γᵢⱼ
  end

  # upper triangular entries
  @inbounds for j in 1:(n * k)
    for i in 1:(j - 1)
      Γ[i, j] = Γ[j, i]
    end
  end

  Γ
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("variogram/gaussian.jl")
include("variogram/spherical.jl")
include("variogram/exponential.jl")
include("variogram/matern.jl")
include("variogram/cubic.jl")
include("variogram/pentaspherical.jl")
include("variogram/sinehole.jl")
include("variogram/circular.jl")
include("variogram/power.jl")
include("variogram/nugget.jl")
