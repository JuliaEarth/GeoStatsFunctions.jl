# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _fit(
  T::Type{<:PiecewiseLinearTransiogram},
  t::EmpiricalTransiogram,
  ::FitAlgo;
)
  T(t.abscissas, t.ordinates)
end
