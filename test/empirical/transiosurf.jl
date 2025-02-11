@testset "EmpiricalTransiogramSurface" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  data = georef((c=[v < 0 ? 1 : 2 for v in img],))
  t = EmpiricalTransiogramSurface(data, :c, maxlag=50.0)
  @test sprint(show, t) == "EmpiricalTransiogramSurface(nangs: 100, nlags: 20)"
  @test sprint(show, MIME"text/plain"(), t) == """
  EmpiricalTransiogramSurface
  ├─ nangs: 100
  └─ nlags: 20"""
end
