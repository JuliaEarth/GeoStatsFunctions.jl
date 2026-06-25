# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

using GeoTables
using TableTransforms

@setup_workload begin
  # synthetic data
  d = georef((z=[sin(i/2) + sin(j/2) for i in 1:50, j in 1:50],))
  c = d |> Indicator("z", k=2) |> Select(1 => "z")

  # empirical variograms and transiograms
  g = EmpiricalVariogram(d, "z", maxlag=25)
  t = EmpiricalTransiogram(c, "z", maxlag=25)

  # theoretical variograms and transiograms
  γ = ExponentialVariogram(ranges=(3, 2, 1))
  τ = ExponentialTransiogram(ranges=(3, 2, 1))

  @compile_workload begin
    for f in [g, t, γ, τ]
      funplot(f)
    end
  end
end
