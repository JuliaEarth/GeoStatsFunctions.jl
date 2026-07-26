@testset "variogramsurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  gtb = georef((z=img,))
  g = variogramsurface(gtb, "z", maxlag=50.0)
  @test sprint(show, g) == "EmpiricalVariogramSurface"
  @test sprint(show, MIME"text/plain"(), g) == """
  EmpiricalVariogramSurface
  ├─ nangs: 50
  └─ nlags: 20"""
end
