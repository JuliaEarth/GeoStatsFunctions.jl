@testset "EmpiricalVarioplane" begin
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  data = georef((z=img,))
  g = EmpiricalVarioplane(data, :z, maxlag=50.0)
  @test sprint(show, g) == "EmpiricalVarioplane"
  @test sprint(show, MIME"text/plain"(), g) == """
  EmpiricalVarioplane
    50 angles
    └─0.00°
    └─3.67°
    └─7.35°
    └─11.02°
    └─14.69°
    ⋮
    └─165.31°
    └─168.98°
    └─172.65°
    └─176.33°
    └─180.00°"""
end