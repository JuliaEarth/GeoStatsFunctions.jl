@testset "Covariance" begin
  x, y = rand(Point), rand(Point)
  for (CovType, VarioType) in [
    (CircularCovariance, CircularVariogram),
    (CubicCovariance, CubicVariogram),
    (ExponentialCovariance, ExponentialVariogram),
    (GaussianCovariance, GaussianVariogram),
    (MaternCovariance, MaternVariogram),
    (PentasphericalCovariance, PentasphericalVariogram),
    (SineHoleCovariance, SineHoleVariogram),
    (SphericalCovariance, SphericalVariogram)
  ]
    γ = VarioType()
    cov = CovType(γ)
    @test cov(x, y) == sill(γ) - γ(x, y)
  end

  for (CovType, VarioType) in [
    (CircularCovariance, CircularVariogram),
    (CubicCovariance, CubicVariogram),
    (ExponentialCovariance, ExponentialVariogram),
    (GaussianCovariance, GaussianVariogram),
    (MaternCovariance, MaternVariogram),
    (PentasphericalCovariance, PentasphericalVariogram),
    (SineHoleCovariance, SineHoleVariogram),
    (SphericalCovariance, SphericalVariogram)
  ]
    γ = VarioType(sill=1.5)
    cov = CovType(γ)
    @test cov(x, y) == 1.5 - γ(x, y)
  end

  𝒟 = PointSet([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
  Γ = GeoStatsFunctions.pairwise(GaussianCovariance(), 𝒟)
  @test eltype(Γ) == Float64
  @test issymmetric(Γ)

  𝒟 = PointSet([(1.0f0, 0.0f0, 0.0f0), (0.0f0, 1.0f0, 0.0f0), (0.0f0, 0.0f0, 1.0f0)])
  Γ_f = GeoStatsFunctions.pairwise(GaussianCovariance(range=1.0f0, sill=1.0f0, nugget=0.0f0), 𝒟)
  @test eltype(Γ_f) == Float32
  @test issymmetric(Γ_f)

  cov = GaussianCovariance(range=1.0, sill=1.0, nugget=0.0)
  @test isstationary(cov)
  @test isisotropic(cov)
  @test sill(cov) == 1.0
  @test nugget(cov) == 0.0
  @test metricball(cov) == MetricBall(1.0)
  @test range(cov) == 1.0u"m"
  @test GeoStatsFunctions.scale(cov, 2) == GaussianCovariance(range=2.0, sill=1.0, nugget=0.0)

  # shows
  cov = CircularCovariance()
  @test sprint(show, cov) == "CircularCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CircularCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = CubicCovariance()
  @test sprint(show, cov) == "CubicCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CubicCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = ExponentialCovariance()
  @test sprint(show, cov) == "ExponentialCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  ExponentialCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = GaussianCovariance()
  @test sprint(show, cov) == "GaussianCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  GaussianCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = MaternCovariance()
  @test sprint(show, cov) == "MaternCovariance(sill: 1.0, nugget: 0.0, order: 1.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  MaternCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ order: 1.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = PentasphericalCovariance()
  @test sprint(show, cov) == "PentasphericalCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  PentasphericalCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = SineHoleCovariance()
  @test sprint(show, cov) == "SineHoleCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SineHoleCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""

  cov = SphericalCovariance()
  @test sprint(show, cov) == "SphericalCovariance(sill: 1.0, nugget: 0.0, range: 1.0 m, distance: Euclidean)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SphericalCovariance
  ├─ sill: 1.0
  ├─ nugget: 0.0
  ├─ range: 1.0 m
  └─ distance: Euclidean"""
end
