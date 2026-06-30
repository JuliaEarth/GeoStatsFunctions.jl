# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

aslen(x::Len) = x
aslen(x::Number) = x * u"m"
aslen(::Quantity) = throw(ArgumentError("invalid length unit"))

_maxlag(f::GeoStatsFunction) = 3range(f)
_maxlag(::PowerVariogram) = 3.0u"m"
_maxlag(::NuggetEffect) = 3.0u"m"
_maxlag(f::EmpiricalGeoStatsFunction) = last(f.abscissas)
_maxlag(f::EmpiricalGeoStatsSurface) = last(f.rs)

_ylabel(f::GeoStatsFunction) = "function value"
_ylabel(f::Variogram) = "variogram"
_ylabel(f::Covariance) = "covariance"
_ylabel(f::Transiogram) = "probability"
_ylabel(f::CompositeFunction) = _ylabel(first(last(structures(f))))
_ylabel(f::EmpiricalGeoStatsFunction) = "function value"
_ylabel(f::EmpiricalVariogram) = "variogram"
_ylabel(f::EmpiricalTransiogram) = "probability"

_eval(f, hs) = isisotropic(f) ? _isoeval(f, hs) : _anisoeval(f, hs)

function _isoeval(f, hs)
  # auxiliary parameters
  n = nvariables(f)

  # evaluate all lags
  fs = f.(hs)

  # reshape result
  Fᵢₛₒ = [getindex.(fs, i, j) for i in 1:n, j in 1:n]

  [Fᵢₛₒ]
end

function _anisoeval(f, hs)
  # auxiliary parameters
  n = nvariables(f)

  # reference point and basis vectors
  p, v = _anisobasis(f)

  # evaluate along basis vectors
  map(v) do vⱼ
    fs = [f(p, p + ustrip(h) * vⱼ) for h in hs]
    [getindex.(fs, i, j) for i in 1:n, j in 1:n]
  end
end

function _anisobasis(f)
  # auxiliary parameters
  b = metricball(f)
  R = rotation(b)
  r = radii(b)
  d = length(r)
  U = eltype(r)

  # reference point
  p = Point(ntuple(i -> U(0), d))

  # rotated basis
  v = ntuple(d) do j
    R * Vec(ntuple(i -> U(i == j), d))
  end

  p, v
end

function _anisobasis(f::CarleTransiogram{N}) where {N}
  # auxiliary parameters
  U = typeof(range(f))

  # reference point
  p = Point(ntuple(i -> U(0), N))

  # axis-aligned basis
  v = ntuple(N) do j
    Vec(ntuple(i -> U(i == j), N))
  end

  p, v
end
