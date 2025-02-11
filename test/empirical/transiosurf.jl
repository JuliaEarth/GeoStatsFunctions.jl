@testset "EmpiricalTransiogramSurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  data = georef((c=[v < 0 ? 1 : 2 for v in img],))
  t = EmpiricalTransiogramSurface(data, :c, maxlag=50.0)
  @test sprint(show, t) == "EmpiricalTransiogramSurface"
  @test sprint(show, MIME"text/plain"(), t) == """
  EmpiricalTransiogramSurface
    100 angles
    └─0.00°
    └─3.64°
    └─7.27°
    └─10.91°
    └─14.55°
    ⋮
    └─345.45°
    └─349.09°
    └─352.73°
    └─356.36°
    └─360.00°"""
end
