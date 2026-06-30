# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    funplot(f; [options])

Plot the geostatistical function `f` with given `options`.

## Common options:

* `color`     - color of function graph
* `linewidth` - line width of function graph
* `maxlag`    - maximum lag distance

## Empirical function options:

* `linestyle` - line style of function graph
* `pointsize` - size of points
* `showtext`  - show text counts
* `textsize`  - size of text counts
* `showhist`  - show histogram
* `histcolor` - color of histogram

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

"""
    hscatter(geotable, vars; [options])

h-scatter plot for all variables `vars` stored in the `geotable`.
Optionally, specify the options documented below.

## Options

* `lag`      - lag distance in length units (default to `0.0u"m"`)
* `tol`      - tolerance for lag distance (default to `0.1u"m"`)
* `distance` - distance from Distances.jl (default to `Euclidean()`)
* `nmax`     - maximum number of samples (default to `4000`)
* `size`     - size of points (default to `2`)
* `color`    - color of points (default to `:black`)
* `alpha`    - transparency of points (default to `1.0`)
* `rcolor`   - color of regression line (default to `:salmon`)
* `icolor`   - color of identity line (default to `:black`)
* `ccolor`   - color of center lines (default to `:teal`)

## Examples

```julia
# h-scatter of z vs. z at lag 1.0
hscatter(geotable, "z", lag=1.0)

# h-scatter of z vs. w at lag 2.0
hscatter(geotable, ["z", "w"], lag=2.0)
```

### Notes

This function will only work in the presence of
a Makie.jl backend via package extensions in
Julia v1.9 or later versions of the language.
"""
function hscatter end
