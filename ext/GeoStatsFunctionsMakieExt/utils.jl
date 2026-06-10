# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const Len{T} = Quantity{T,u"𝐋"}

_maxlag(f::GeoStatsFunction) = 3range(f)
_maxlag(::PowerVariogram) = 3.0u"m"
_maxlag(::NuggetEffect) = 3.0u"m"
_maxlag(f::EmpiricalGeoStatsFunction) = last(f.abscissas)
_maxlag(f::EmpiricalGeoStatsSurface) = last(f.rs)

_istransiogram(f) = false
_istransiogram(t::Transiogram) = true
_istransiogram(t::EmpiricalTransiogram) = true
_istransiogram(t::EmpiricalTransiogramSurface) = true

_layout(fig::Makie.Figure) = fig.layout
_layout(gl::Makie.GridLayout) = gl
_layout(gp::Makie.GridPosition) = Makie.GridLayout(gp)
_layout(gp::Makie.GridSubposition) = Makie.GridLayout(gp)