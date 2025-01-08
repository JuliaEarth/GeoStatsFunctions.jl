@testset "Sampling" begin
  γ = GaussianVariogram()
  seg = Segment((0.0, 0.0), (1.0, 1.0))
  ps = GeoStatsFunctions._sample(γ, seg) |> collect
  @test all(p -> Point(0.0, 0.0) ≤ p ≤ Point(1.0, 1.0), ps)
  @test length(ps) == 3

  γ = GaussianVariogram()
  quad = Quadrangle((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
  ps = GeoStatsFunctions._sample(γ, quad) |> collect
  @test all(p -> Point(0.0, 0.0) ≤ p ≤ Point(1.0, 1.0), ps)
  @test length(ps) == 3 * 3

  γ = GaussianVariogram()
  hex = Hexahedron(
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0)
  )
  ps = GeoStatsFunctions._sample(γ, hex) |> collect
  @test all(p -> Point(0.0, 0.0, 0.0) ≤ p ≤ Point(1.0, 1.0, 1.0), ps)
  @test length(ps) == 3 * 3 * 3

  # deterministic samples in arbitrary geometries
  γ = GaussianVariogram()
  G = PolyArea((0.0, 0.0), (0.5, -1.5), (1.0, 0.0), (1.5, 0.5), (1.0, 1.0), (0.5, 1.5), (-0.5, 0.5))
  ps1 = GeoStatsFunctions._sample(γ, G)
  ps2 = GeoStatsFunctions._sample(γ, G)
  @test ps1 == ps2

  # samples with nugget effect model
  γ = NuggetEffect()
  h = Hexahedron(
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0)
  )
  ps = GeoStatsFunctions._sample(γ, h) |> collect
  @test length(ps) == 3 * 3 * 3

  # minimum allocation while sampling points
  γ = GaussianVariogram()
  p = Point(0.0, 0.0, 0.0)
  ps = GeoStatsFunctions._sample(γ, p) |> collect
  @test ps == [p]
  @test (@allocated GeoStatsFunctions._sample(γ, p)) ≤ 32
end
