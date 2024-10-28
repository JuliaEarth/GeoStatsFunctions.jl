# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _fit(T::Type{<:PiecewiseLinearTransiogram}, t::EmpiricalTransiogram, ::FitAlgo)
  τ = T(t.abscissas, t.ordinates)
  ϵ = 0.0
  τ, ϵ
end
