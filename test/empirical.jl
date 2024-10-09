@testset "Empirical" begin
  @testset "Variogram" begin
    # homogeneous field has zero variogram
    sdata = georef((z=ones(3),), [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
    γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2.0)
    x, y, n = values(γ)
    @test x ≈ [1 / 2, √2] * u"m"
    @test y[2] == 0.0
    @test n == [0, 3]

    # basic test on number of lags
    sdata = georef((z=[1.0, 0.0, 1.0],), [(25.0, 25.0), (50.0, 75.0), (75.0, 50.0)])
    γ = EmpiricalVariogram(sdata, :z, nlags=20, maxlag=1.0)
    x, y, n = values(γ)
    @test length(x) == 20
    @test length(y) == 20
    @test length(n) == 20

    # empirical variogram on integer coordinates
    sdata = georef((z=ones(3),), [(1, 0, 0), (0, 1, 0), (0, 0, 1)])
    γ = EmpiricalVariogram(sdata, :z, nlags=2, maxlag=2, algorithm=:full)
    x, y, n = values(γ)
    @test x ≈ [1 / 2, √2] * u"m"
    @test y[2] == 0.0
    @test n == [0, 3]

    # empirical variogram with only missing data
    z = Union{Float64,Missing}[missing, missing, missing]
    𝒟 = georef((z=z,), rand(Point, 3))
    γ = EmpiricalVariogram(𝒟, :z, maxlag=1.0, nlags=5)
    x, y, n = values(γ)
    @test x == [0.1, 0.3, 0.5, 0.7, 0.9] * u"m"
    @test all(iszero.(n))

    # accumulation algorithms give the same result
    rng = StableRNG(123)
    sdata = georef((z=rand(rng, 1000),), rand(rng, Point, 1000))
    γ₁ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algorithm=:full)
    γ₂ = EmpiricalVariogram(sdata, :z, maxlag=0.01, algorithm=:ball)
    @test isequal(values(γ₁), values(γ₂))

    # custom distance is recorded
    rng = StableRNG(123)
    sdata = georef((z=rand(rng, 2),), [Point(LatLon(0.0, 0.0)), Point(LatLon(0.0, 90.0))])
    γ = EmpiricalVariogram(sdata, :z, distance=Haversine(6371.0), algorithm=:full)
    @test distance(γ) == Haversine(6371.0)

    # print methods
    rng = StableRNG(123)
    d = georef((z=rand(rng, 100, 100),))
    γ = EmpiricalVariogram(d, :z)
    @test sprint(show, γ) ==
          "EmpiricalVariogram(abscissa: [0.25 m, ..., 9.93304 m], ordinate: [0.0, ..., 0.0841979], distance: Euclidean(0.0), estimator: MatheronEstimator(), npairs: 1447200)"
    @test sprint(show, MIME"text/plain"(), γ) == """
    EmpiricalVariogram
    ├─ abscissa: [0.25 m, 1.0 m, 1.41421 m, ..., 8.7407 m, 9.28182 m, 9.93304 m]
    ├─ ordinate: [0.0, 0.0843099, 0.0845995, ..., 0.0838336, 0.0839823, 0.0841979]
    ├─ distance: Euclidean(0.0)
    ├─ estimator: MatheronEstimator()
    └─ npairs: 1447200"""

    # test variography with compositional data
    data = georef((z=rand(Composition{3}, 100),), rand(Point, 100))
    γ = EmpiricalVariogram(data, :z, maxlag=1.0, algorithm=:full)
    x, y, n = values(γ)
    @test all(≥(0u"m"), x)
    @test all(≥(0), y)
    @test all(≥(0), n)

    # test variography with unitful data
    data = georef((z=[1 * u"K" for i in 1:100],), rand(Point, 100))
    γ = EmpiricalVariogram(data, :z, nlags=20)
    x, y, n = values(γ)
    @test all(≥(0u"m"), x)
    @test y == fill(0.0 * u"K^2", 20)

    # Matheron's vs Cressie's estimator
    img = readdlm(joinpath(datadir, "Gaussian30x10.txt"))
    data = georef((; Z=img))
    γ₁ = EmpiricalVariogram(data, :Z, maxlag=50.0, estimator=:matheron)
    γ₂ = EmpiricalVariogram(data, :Z, maxlag=50.0, estimator=:cressie)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test all(isapprox.(y₁, y₂, atol=0.1))
    @test n₁ == n₂

    # specify variables as strings
    img = readdlm(joinpath(datadir, "Gaussian30x10.txt"))
    data = georef((; Z=img))
    γ = EmpiricalVariogram(data, "Z", maxlag=50.0)
    x, y, n = values(γ)
    @test all(≥(0u"m"), x)
    @test all(>(0.8), y[11:end])
    @test all(≥(0), n)
  end

  @testset "Varioplane" begin
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    data = georef((z=img,))
    γ = EmpiricalVarioplane(data, :z, maxlag=50.0)
    @test sprint(show, γ) == "EmpiricalVarioplane"
    @test sprint(show, MIME"text/plain"(), γ) == """
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

  @testset "Directional" begin
    # merge operation does not produce NaN
    dir = (0.286788, -0.496732, -0.819152)
    𝒟 = georef(CSV.File(joinpath(datadir, "nanlags.csv")), (:X, :Y, :Z))
    γ = DirectionalVariogram(dir, 𝒟, :Cu, dtol=45, maxlag=150, nlags=20)
    x, y, n = values(γ)
    @test !any(isnan.(x))
    @test !any(isnan.(y))
    @test !any(isnan.(n))

    # directional variogram and known anisotropy ratio
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    sdata = georef((z=img,))
    γhor = DirectionalVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γver = DirectionalVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γₕ = GeoStatsFunctions.fit(GaussianVariogram, γhor)
    γᵥ = GeoStatsFunctions.fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
  end

  @testset "Planar" begin
    # directional equals planar rotated by 90 degrees in 2D
    img = readdlm(joinpath(datadir, "anisotropic.tsv"))
    sdata = georef((z=img,))
    γ₁ = PlanarVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γ₂ = DirectionalVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂
    γ₁ = PlanarVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γ₂ = DirectionalVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    x₁, y₁, n₁ = values(γ₁)
    x₂, y₂, n₂ = values(γ₂)
    @test x₁ == x₂
    @test y₁ ≈ y₂
    @test n₁ == n₂

    # planar variogram and known anisotropy ratio
    γhor = PlanarVariogram((0.0, 1.0), sdata, :z, maxlag=50.0)
    γver = PlanarVariogram((1.0, 0.0), sdata, :z, maxlag=50.0)
    γₕ = GeoStatsFunctions.fit(GaussianVariogram, γhor)
    γᵥ = GeoStatsFunctions.fit(GaussianVariogram, γver)
    @test range(γₕ) / range(γᵥ) ≈ 3.0 atol = 0.1
  end
end
