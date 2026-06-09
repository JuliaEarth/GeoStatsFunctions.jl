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

const GridPos = Union{Makie.Figure, Makie.GridLayout, Makie.GridPosition, Makie.GridSubposition}

function get_layout(gp::Makie.Figure)
  gp.layout
end

function get_layout(gp::Makie.GridLayout)
  gp
end

function get_layout(gp::Union{Makie.GridPosition, Makie.GridSubposition})
  Makie.GridLayout(gp)
end

