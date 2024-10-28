# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    fit(T, t, algo=WeightedLeastSquares())

Fit theoretical transiogram type `T` to empirical transiogram `t` using algorithm `algo`.

## Examples

```julia
julia> fit(PiecewiseLinearTransiogram, t)
```
"""
fit(::Type{PiecewiseLinearTransiogram}, t::EmpiricalTransiogram) =
  PiecewiseLinearTransiogram(t.abscissas, t.ordinates)
