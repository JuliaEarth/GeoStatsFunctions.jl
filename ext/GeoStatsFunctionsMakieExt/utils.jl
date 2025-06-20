# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

_maxlag(f::GeoStatsFunction) = 3range(f)
_maxlag(::PowerVariogram) = 3.0u"m"
_maxlag(::NuggetEffect) = 3.0u"m"
_maxlag(f::EmpiricalGeoStatsFunction) = last(f.abscissas)
_maxlag(f::EmpiricalGeoStatsSurface) = last(f.rs)

_addunit(x::Number, u) = x * u
_addunit(x::Len, _) = x
_addunit(x::Quantity, _) = throw(ArgumentError("$(unit(x)) is not a valid length unit"))

_istransiogram(f) = false
_istransiogram(t::Transiogram) = true
_istransiogram(t::EmpiricalTransiogram) = true
_istransiogram(t::EmpiricalTransiogramSurface) = true
