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

* `style`       - style (line style)
* `pointsize`   - size of points
* `showtext`    - show text counts
* `textsize`    - size of text counts
* `showhist`    - show histogram
* `histcolor`   - color of histogram

### Notes

This function will only work in the presence of
a Makie.jl backend via package extensions in
Julia v1.9 or later versions of the language.
"""
function funplot end

"""
    funplot!(fig, f; [options])

Mutating version of [`funplot`[@ref] where the figure `fig`
is updated with the plot of the geostatistical function `f`.

## Examples

```julia
# initialize figure with Gaussian variogram
fig = funplot(GaussianVariogram())

# add spherical variogram to figure
funplot!(fig, SphericalVariogram())
```

See the documentation of [`funplot`](@ref) for `options`.
"""
function funplot! end

"""
    surfplot(f; [options])

Plot the geostatistical surface `f` with given `options`.

## Common options

* `colormap` - Color map
* `maxlag`   - maximum lag
* `labels`   - variable names

## Theoretical function options

* `normal` - Normal direction to plane (default to vertical)
* `nlags`  - Number of lags (default to `20`)
* `nangs`  - Number of angles (default to `50`)

## Examples

```julia
# initialize figure with Gaussian variogram
fig = surfplot(GaussianVariogram())

# add spherical variogram to figure
surfplot!(fig, SphericalVariogram())
```

### Notes

This function will only work in the presence of
a Makie.jl backend via package extensions in
Julia v1.9 or later versions of the language.
"""
function surfplot end

"""
    surfplot!(fig, f; [options])

Mutating version of [`surfplot`[@ref] where the figure `fig`
is updated with the plot of the geostatistical surface `f`.

See the documentation of [`surfplot`](@ref) for `options`.
"""
function surfplot! end
