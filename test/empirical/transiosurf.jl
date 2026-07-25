@testset "EmpiricalTransiogramSurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  gtb = georef((c=[v < 0 ? 1 : 2 for v in img],))
  t = transiogramsurface(gtb, "c", maxlag=50.0)
  @test sprint(show, t) == "EmpiricalTransiogramSurface"
  @test sprint(show, MIME"text/plain"(), t) == """
  EmpiricalTransiogramSurface
  ├─ nangs: 100
  └─ nlags: 20"""
end
