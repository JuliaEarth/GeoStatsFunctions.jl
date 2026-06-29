@testset "EmpiricalVariogram" begin
  # homogeneous field has zero variogram
  gtb = georef((z=ones(3),), [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
  g = EmpiricalVariogram(gtb, "z", nlags=2, maxlag=2.0)
  @test g.abscissas ≈ [1 / 2, √2] * u"m"
  @test g.ordinates[1][2] == 0.0
  @test g.counts == [0, 3]

  # basic test on number of lags
  gtb = georef((z=[1.0, 0.0, 1.0],), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
  g = EmpiricalVariogram(gtb, "z", nlags=20, maxlag=1.0)
  @test length(g.abscissas) == 20
  @test length(g.ordinates) == 1
  @test length(g.ordinates[1]) == 20
  @test length(g.counts) == 20

  # empirical variogram on integer coordinates
  gtb = georef((z=ones(3),), [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
  g = EmpiricalVariogram(gtb, "z", nlags=2, maxlag=2)
  @test g.abscissas ≈ [1 / 2, √2] * u"m"
  @test g.ordinates[1][2] == 0.0
  @test g.counts == [0, 3]

  # empirical variogram with only missing data
  z = Union{Float64,Missing}[missing, missing, missing]
  gtb = georef((z=z,), [(0.0, 0.0), (0.1, 0.0), (0.5, 0.0)])
  g = EmpiricalVariogram(gtb, "z", maxlag=1.0, nlags=5)
  @test g.abscissas == [0.1, 0.4, 0.5, 0.7, 0.9] * u"m"
  @test all(ismissing, g.ordinates[1])
  @test g.counts == [1, 1, 1, 0, 0]

  # lag search methods give the same result
  rng = StableRNG(123)
  gtb = georef((z=rand(rng, 1000),), rand(rng, Point, 1000))
  g₁ = EmpiricalVariogram(gtb, "z", maxlag=0.01, lagsearch=:full)
  g₂ = EmpiricalVariogram(gtb, "z", maxlag=0.01, lagsearch=:ball)
  @test isequal(g₁.abscissas, g₂.abscissas)
  @test isequal(g₁.ordinates, g₂.ordinates)
  @test isequal(g₁.counts, g₂.counts)

  # custom distance is recorded
  rng = StableRNG(123)
  gtb = georef((z=rand(rng, 2),), [Point(LatLon(0.0, 0.0)), Point(LatLon(0.0, 90.0))])
  g = EmpiricalVariogram(gtb, "z", distance=Haversine(6371.0), lagsearch=:full)
  @test g.distance == Haversine(6371.0)

  # test variography with compositional data
  gtb = georef((z=rand(Composition{3}, 100),), rand(Point, 100))
  g = EmpiricalVariogram(gtb, "z", maxlag=1.0)
  @test all(≥(0u"m"), g.abscissas)
  @test all(≥(0), g.ordinates[1])
  @test all(≥(0), g.counts)

  # test variography with unitful data
  gtb = georef((z=[1 * u"K" for i in 1:100],), rand(Point, 100))
  g = EmpiricalVariogram(gtb, "z", nlags=20)
  @test all(≥(0u"m"), g.abscissas)
  @test g.ordinates[1] == fill(0.0 * u"K^2", 20)

  # Matheron's vs Cressie's estimator
  img = readdlm(joinpath(datadir, "Gaussian30x10.txt"))
  gtb = georef((; z=img))
  g₁ = EmpiricalVariogram(gtb, "z", maxlag=50.0, estimator=:matheron)
  g₂ = EmpiricalVariogram(gtb, "z", maxlag=50.0, estimator=:cressie)
  @test g₁.abscissas == g₂.abscissas
  @test all(isapprox.(g₁.ordinates[1], g₂.ordinates[1], atol=0.1))
  @test g₁.counts == g₂.counts

  # specify variables as strings
  img = readdlm(joinpath(datadir, "Gaussian30x10.txt"))
  gtb = georef((; z=img))
  g = EmpiricalVariogram(gtb, "z", maxlag=50.0)
  @test all(≥(0u"m"), g.abscissas)
  @test all(>(0.8), g.ordinates[1][11:end])
  @test all(≥(0), g.counts)
end

@testset "DirectionalVariogram" begin
  # merge operation does not produce NaN
  dir = (0.286788, -0.496732, -0.819152)
  gtb = georef(CSV.File(joinpath(datadir, "nanlags.csv")), (:X, :Y, :Z))
  g = DirectionalVariogram(dir, gtb, "Cu", dtol=45, maxlag=150, nlags=20)
  @test !any(isnan.(g.abscissas))
  @test !any(isnan.(g.ordinates[1]))
  @test !any(isnan.(g.counts))

  # directional variogram and known anisotropy ratio
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  gtb = georef((z=img,))
  gₕ = DirectionalVariogram((1.0, 0.0), gtb, "z", maxlag=50.0)
  gᵥ = DirectionalVariogram((0.0, 1.0), gtb, "z", maxlag=50.0)
  γₕ = GeoStatsFunctions.fit(GaussianVariogram, gₕ)
  γᵥ = GeoStatsFunctions.fit(GaussianVariogram, gᵥ)
  @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
end

@testset "PlanarVariogram" begin
  # directional equals planar rotated by 90 degrees in 2D
  img = readdlm(joinpath(datadir, "anisotropic.tsv"))
  gtb = georef((z=img,))
  g₁ = PlanarVariogram((0.0, 1.0), gtb, "z", maxlag=50.0)
  g₂ = DirectionalVariogram((1.0, 0.0), gtb, "z", maxlag=50.0)
  @test g₁.abscissas == g₂.abscissas
  @test g₁.ordinates[1] ≈ g₂.ordinates[1]
  @test g₁.counts == g₂.counts
  g₁ = PlanarVariogram((1.0, 0.0), gtb, "z", maxlag=50.0)
  g₂ = DirectionalVariogram((0.0, 1.0), gtb, "z", maxlag=50.0)
  @test g₁.abscissas == g₂.abscissas
  @test g₁.ordinates[1] ≈ g₂.ordinates[1]
  @test g₁.counts == g₂.counts

  # planar variogram and known anisotropy ratio
  gₕ = PlanarVariogram((0.0, 1.0), gtb, "z", maxlag=50.0)
  gᵥ = PlanarVariogram((1.0, 0.0), gtb, "z", maxlag=50.0)
  γₕ = GeoStatsFunctions.fit(GaussianVariogram, gₕ)
  γᵥ = GeoStatsFunctions.fit(GaussianVariogram, gᵥ)
  @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
end
