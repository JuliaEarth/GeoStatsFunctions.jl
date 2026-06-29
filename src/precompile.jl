# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  gtb2D = georef((; z=rand(5, 5)))
  gtb3D = georef((; z=rand(5, 5, 5)))
  @compile_workload begin
    g = EmpiricalVariogram(gtb2D, "z")
    g = EmpiricalVariogram(gtb3D, "z")
    γ = GeoStatsFunctions.fit(Variogram, g)
  end
end
