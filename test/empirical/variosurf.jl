@testset "variogramsurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  gtb = georef((z=img,))
  g = variogramsurface(gtb, "z", maxlag=50.0)
  @test sprint(show, g) == "EmpiricalVariogramSurface"
  @test sprint(show, MIME"text/plain"(), g) == """
  EmpiricalVariogramSurface
  ├─ nangs: 50
  └─ nlags: 20"""
  @test length(g.ns) == 50 # count for each of default nangs
  @test all([length(n) == length(g.rs) for n in g.ns]) # count for every abscissa
end
