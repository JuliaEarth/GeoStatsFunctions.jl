# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    varioplot(γ; [options])

Plot the variogram or varioplane `γ` with given `options`.

## Empirical variogram options:

* `color`       - color of variogram
* `pointsize`   - size of points of variogram
* `segmentsize` - size of segments of variogram
* `showtext`    - show text counts
* `textsize`    - size of text counts
* `showhist`    - show histogram
* `histcolor`   - color of histogram

## Empirical varioplane options:

* `colormap`   - color map of varioplane
* `showrange`  - show varioplane range
* `rangecolor` - color of varioplane range
* `rangemodel` - range model (e.g. `SphericalVariogram`)

## Theoretical variogram options:

* `maxlag` - maximum lag for theoretical model

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function varioplot end
function varioplot! end
