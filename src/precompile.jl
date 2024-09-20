# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

using PrecompileTools

@setup_workload begin
  pts3D = rand(Point, 100)
  pts2D = rand(Point, 100) |> Shadow("xy")
  gtb3D = georef((; Z=rand(100)), pts3D)
  gtb2D = georef((; Z=rand(100)), pts2D)
  @compile_workload begin
    g = EmpiricalVariogram(gtb3D, "Z")
    g = EmpiricalVariogram(gtb2D, "Z")
    γ = GeoStatsFunctions.fit(Variogram, g)
  end
end
