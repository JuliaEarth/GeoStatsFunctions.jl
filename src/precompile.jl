# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  t1 = georef((; z=rand(5, 5)))
  t2 = georef((; z=rand(5, 5, 5)))
  Gs = [GaussianVariogram, SphericalVariogram, ExponentialVariogram]
  @compile_workload begin
    g1 = EmpiricalVariogram(t1, "z")
    g2 = EmpiricalVariogram(t2, "z")
    γ1 = GeoStatsFunctions.fit(Gs, g1)
    γ2 = GeoStatsFunctions.fit(Gs, g2)
  end
end
