# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    varioplot(γ; [options])

Plot the variogram or varioplane `γ` with given `options`.

## Common variogram options:

* `color`  - color of variogram
* `size`   - size of variogram
* `maxlag` - maximum lag of variogram

## Empirical variogram options:

* `pointsize`   - size of points of variogram
* `showtext`    - show text counts
* `textsize`    - size of text counts
* `showhist`    - show histogram
* `histcolor`   - color of histogram

## Empirical varioplane options:

* `colormap` - color map of varioplane

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function varioplot end
function varioplot! end

"""
    transioplot(t; [options])

Plot the transiogram `t` with given `options`.

## Common transiogram options:

* `color`  - color of transiogram
* `size`   - size of transiogram
* `maxlag` - maximum lag of variogram
* `levels` - categorical levels

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function transioplot end
function transioplot! end
