@testset "EmpiricalVariogramSurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  data = georef((z=img,))
  g = EmpiricalVariogramSurface(data, :z, maxlag=50.0)
  @test sprint(show, g) == "EmpiricalVariogramSurface(nangs: 50, nlags: 20)"
  @test sprint(show, MIME"text/plain"(), g) == """
  EmpiricalVariogramSurface
  ├─ nangs: 50
  └─ nlags: 20"""
end
