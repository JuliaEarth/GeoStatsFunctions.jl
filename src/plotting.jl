# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    funplot(f; [options])

Plot the geostatistical function `f` with given `options`.

## Common options:

* `color`  - color
* `size`   - size (line width)
* `maxlag` - maximum lag
* `labels` - variable names

## Empirical function options:

* `pointsize`   - size of points
* `showtext`    - show text counts
* `textsize`    - size of text counts
* `showhist`    - show histogram
* `histcolor`   - color of histogram

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function funplot end

"""
    planeplot(f; [options])

Plot the varioplane or transioplane `f` with given `options`.

## Common options

* `colormap` - Color map
* `maxlag`   - maximum lag
* `labels`   - variable names

## Theoretical function options

* `normal` - Normal direction to plane (default to vertical)
* `nlags`  - Number of lags (default to `20`)
* `nangs`  - Number of angles (default to `50`)

### Notes

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function planeplot end
