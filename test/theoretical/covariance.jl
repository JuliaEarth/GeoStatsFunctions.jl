@testset "Covariance" begin
  x, y = rand(Point), rand(Point)
  for (CovType, VarioType) in [
    (CircularCovariance, CircularVariogram),
    (CubicCovariance, CubicVariogram),
    (ExponentialCovariance, ExponentialVariogram),
    (GaussianCovariance, GaussianVariogram),
    (MaternCovariance, MaternVariogram),
    (PentaSphericalCovariance, PentaSphericalVariogram),
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
    (PentaSphericalCovariance, PentaSphericalVariogram),
    (SineHoleCovariance, SineHoleVariogram),
    (SphericalCovariance, SphericalVariogram)
  ]
    γ = VarioType(sill=1.5)
    cov = CovType(γ)
    @test cov(x, y) == 1.5 - γ(x, y)
  end

  𝒟 = PointSet([(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)])
  C = GeoStatsFunctions.pairwise(GaussianCovariance(), 𝒟)
  @test eltype(C) == Float64
  @test issymmetric(C)

  𝒟 = PointSet([(1.0f0, 0.0f0, 0.0f0), (0.0f0, 1.0f0, 0.0f0), (0.0f0, 0.0f0, 1.0f0)])
  C_f = GeoStatsFunctions.pairwise(GaussianCovariance(range=1.0f0, sill=1.0f0, nugget=0.0f0), 𝒟)
  @test eltype(C_f) == Float32
  @test issymmetric(C_f)

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
  @test sprint(show, cov) == "CircularCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CircularCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = CubicCovariance()
  @test sprint(show, cov) == "CubicCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  CubicCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = ExponentialCovariance()
  @test sprint(show, cov) == "ExponentialCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  ExponentialCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = GaussianCovariance()
  @test sprint(show, cov) == "GaussianCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  GaussianCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = MaternCovariance()
  @test sprint(show, cov) == "MaternCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0, order: 1.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  MaternCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  ├─ nugget: 0.0
  └─ order: 1.0"""

  cov = PentaSphericalCovariance()
  @test sprint(show, cov) == "PentaSphericalCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  PentaSphericalCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = SineHoleCovariance()
  @test sprint(show, cov) == "SineHoleCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SineHoleCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""

  cov = SphericalCovariance()
  @test sprint(show, cov) == "SphericalCovariance(range: 1.0 m, sill: 1.0, nugget: 0.0)"
  @test sprint(show, MIME"text/plain"(), cov) == """
  SphericalCovariance
  ├─ range: 1.0 m
  ├─ sill: 1.0
  └─ nugget: 0.0"""
end
