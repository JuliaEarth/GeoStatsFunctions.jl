@testset "Variogram" begin
  rng = StableRNG(123)
  h = range(0, stop=10, length=50)
  x, y = rand(rng, Point), rand(rng, Point)

  # stationary variogram models
  γs = [
    NuggetEffect(),
    GaussianVariogram(),
    ExponentialVariogram(),
    MaternVariogram(),
    SphericalVariogram(),
    SphericalVariogram(range=2.0),
    CubicVariogram(),
    PentasphericalVariogram(),
    SineHoleVariogram(),
    CircularVariogram()
  ]

  # non-stationary variogram models
  γn = [PowerVariogram(), PowerVariogram(exponent=0.4)]

  # non-decreasing variogram models
  γnd = [
    GaussianVariogram(),
    ExponentialVariogram(),
    MaternVariogram(),
    SphericalVariogram(),
    SphericalVariogram(range=2.0),
    CubicVariogram(),
    PentasphericalVariogram(),
    PowerVariogram()
  ]

  # anisotropic variogram models
  γa = [GaussianVariogram(MetricBall((2.0, 1.0))), MaternVariogram(MetricBall((3.0, 2.0, 1.0)))]

  # check stationarity
  @test all(isstationary, γs)
  @test all(!isstationary, γn)

  # check anisotropy
  @test all(isisotropic, γs)
  @test all(isisotropic, γn)
  @test all(isisotropic, γnd)
  @test isisotropic(sum(γs) + sum(γn) + sum(γnd))
  @test all(!isisotropic, γa)

  # check metric ball
  @test metricball(γa[1]) == MetricBall((2.0, 1.0))
  @test metricball(γa[2]) == MetricBall((3.0, 2.0, 1.0))
  @test metricball(GaussianVariogram(range=2.0)) == MetricBall(2.0)
  @test metricball(SphericalVariogram(range=2.0)) == MetricBall(2.0)
  @test metricball(ExponentialVariogram(range=2.0)) == MetricBall(2.0)
  @test metricball(MaternVariogram(range=2.0)) == MetricBall(2.0)

  # variograms are symmetric under Euclidean distance
  for γ in (γs ∪ γn ∪ γnd ∪ [sum(γs) + sum(γn) + sum(γnd)])
    @test γ(x, y) ≈ γ(y, x)
  end

  # some variograms are non-decreasing
  for γ in (γnd ∪ [sum(γnd)])
    @test all(γ.(h) .≤ γ.(h .+ 1))
  end

  # variograms are valid at the origin
  for γ in (γs ∪ γn ∪ γnd)
    @test !isnan(γ(0.0)) && !isinf(γ(0.0))
  end

  # practical ranges
  for γ in γs
    if !(γ isa NuggetEffect)
      @test isapprox(γ(range(γ)), sill(γ), atol=0.05)
    end
  end

  # nugget effect
  γ = NuggetEffect(nugget=0.2)
  @test nugget(γ) == 0.2
  @test sill(γ) == 0.2
  @test range(γ) == 0u"m"

  # ill-conditioned models and nugget regularization
  # see https://github.com/JuliaEarth/GeoStats.jl/issues/29
  cmat = [
    93.0 90.0 89.0 94.0 93.0 97.0 95.0 88.0 96.0 98.0
    40.0 33.0 34.0 36.0 30.0 39.0 39.0 28.0 25.0 35.0
  ]
  pset = PointSet(Tuple.(eachcol(cmat)))

  # ill-conditioned covariance
  γ = GaussianVariogram(range=20.0)
  C = sill(γ) .- GeoStatsFunctions.pairwise(γ, pset)
  @test cond(C) > 1000.0

  # nugget regularization
  γ = GaussianVariogram(range=20.0, nugget=0.1)
  C = sill(γ) .- GeoStatsFunctions.pairwise(γ, pset)
  @test γ(0) == 0
  @test γ(1e-6) > 0
  @test cond(C) < 100.0

  # sill and nugget in single precision
  for G in [GaussianVariogram, SphericalVariogram, ExponentialVariogram, MaternVariogram]
    γ = G(sill=1.0f0)
    @test typeof(range(γ)) == typeof(1.0u"m")
    @test typeof(sill(γ)) == Float32
    @test typeof(nugget(γ)) == Float32
  end

  # unitful stationary types
  γs = [
    NuggetEffect(1.0u"K^2"),
    GaussianVariogram(sill=1.0u"K^2"),
    ExponentialVariogram(sill=1.0u"K^2"),
    MaternVariogram(sill=1.0u"K^2"),
    SphericalVariogram(sill=1.0u"K^2"),
    CubicVariogram(sill=1.0u"K^2"),
    PentasphericalVariogram(sill=1.0u"K^2"),
    SineHoleVariogram(sill=1.0u"K^2"),
    CircularVariogram(sill=1.0u"K^2")
  ]

  # unitful non-stationary types
  γn = [PowerVariogram(scaling=1.0u"K^2")]
  for γ in γs
    @test unit(γ(1.0)) == u"K^2"
  end

  𝒟 = PointSet([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
  Γ = GeoStatsFunctions.pairwise(GaussianVariogram(), 𝒟)
  @test eltype(Γ) == Float64
  @test issymmetric(Γ)

  𝒟 = PointSet([(1.0f0, 0.0f0, 0.0f0), (0.0f0, 1.0f0, 0.0f0), (0.0f0, 0.0f0, 1.0f0)])
  Γ_f = GeoStatsFunctions.pairwise(GaussianVariogram(range=1.0f0, sill=1.0f0, nugget=0.0f0), 𝒟)
  @test eltype(Γ_f) == Float32
  @test issymmetric(Γ_f)

  𝒟 = CartesianGrid(10, 10)
  Γ = GeoStatsFunctions.pairwise(GaussianVariogram(), view(𝒟, 1:5))
  @test size(Γ) == (5, 5)
  @test issymmetric(Γ)
  Γ = GeoStatsFunctions.pairwise(GaussianVariogram(), view(𝒟, 1:3), view(𝒟, 7:10))
  @test size(Γ) == (3, 4)
  @test all(Γ .> 0)

  # arbitrary collections
  𝒟 = CartesianGrid(10, 10)
  𝒫 = centroid.(𝒟)
  Γ = GeoStatsFunctions.pairwise(GaussianVariogram(), 𝒫)
  @test size(Γ) == (100, 100)
  @test issymmetric(Γ)
  Γ = GeoStatsFunctions.pairwise(GaussianVariogram(), view(𝒫, 1:3), view(𝒫, 7:10))
  @test size(Γ) == (3, 4)
  @test all(Γ .> 0)

  # constructor
  for γ in [
    CircularVariogram(),
    CubicVariogram(),
    ExponentialVariogram(),
    GaussianVariogram(),
    MaternVariogram(),
    NuggetEffect(),
    PentasphericalVariogram(),
    PowerVariogram(),
    SineHoleVariogram(),
    SphericalVariogram()
  ]
    @test GeoStatsFunctions.constructor(γ)() == γ
  end

  # scale
  for γ in [GaussianVariogram(), SphericalVariogram(), ExponentialVariogram()]
    g = GeoStatsFunctions.scale(γ, 2)
    @test range(g) == 2.0u"m"
  end

  # scale with NestedVariogram
  γ = GaussianVariogram(range=2.0) + ExponentialVariogram(range=3.0)
  g = GeoStatsFunctions.scale(γ, 2)
  @test range(g) == 6.0u"m"

  # scale doesn't affect NuggetEffect
  γ = NuggetEffect()
  g = GeoStatsFunctions.scale(γ, 2)
  @test g == γ

  # shows
  γ = CircularVariogram()
  @test sprint(show, γ) == "CircularVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  CircularVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = CubicVariogram()
  @test sprint(show, γ) == "CubicVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  CubicVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = ExponentialVariogram()
  @test sprint(show, γ) == "ExponentialVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  ExponentialVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = GaussianVariogram()
  @test sprint(show, γ) == "GaussianVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  GaussianVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = MaternVariogram()
  @test sprint(show, γ) == "MaternVariogram(sill: 1.0, nugget: 0.0, order: 1.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  MaternVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ order: 1.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = NuggetEffect()
  @test sprint(show, γ) == "NuggetEffect(nugget: 1.0)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  NuggetEffect
  └─ nugget: 1.0"""

  γ = PentasphericalVariogram()
  @test sprint(show, γ) == "PentasphericalVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  PentasphericalVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = PowerVariogram()
  @test sprint(show, γ) == "PowerVariogram(scaling: 1.0, nugget: 0.0, exponent: 1.0)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  PowerVariogram
  ├─ scaling: 1.0
  ├─ nugget: 0.0
  └─ exponent: 1.0"""

  γ = SineHoleVariogram()
  @test sprint(show, γ) == "SineHoleVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  SineHoleVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  γ = SphericalVariogram()
  @test sprint(show, γ) == "SphericalVariogram(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), γ) == """
  SphericalVariogram
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""
end
