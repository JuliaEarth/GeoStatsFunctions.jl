@testset "Fitting variograms" begin
  img = readdlm(joinpath(datadir, "WalkerLake.txt"))
  d = georef((; Z=img))
  g = EmpiricalVariogram(d, :Z, maxlag=15.0)

  # all fits lead to similar sill
  γ₁ = GeoStatsFunctions.fit(GaussianVariogram, g)
  γ₂ = GeoStatsFunctions.fit(SphericalVariogram, g)
  γ₃ = GeoStatsFunctions.fit(ExponentialVariogram, g)
  γ₄ = GeoStatsFunctions.fit(MaternVariogram, g)
  @test isapprox(sill(γ₁), 0.054, atol=1e-3)
  @test isapprox(sill(γ₂), 0.054, atol=1e-3)
  @test isapprox(sill(γ₃), 0.054, atol=1e-3)
  @test isapprox(sill(γ₄), 0.054, atol=1e-3)

  # fix parameters
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, range=12.0u"m")
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, sill=0.07)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, nugget=0.05)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, range=12.0u"m", sill=0.07)
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, range=12.0u"m", nugget=0.05)
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, sill=0.07, nugget=0.05)
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, range=12.0u"m", sill=0.07, nugget=0.05)
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")
  @test isapprox(sill(γ), 0.07, atol=1e-3)
  @test isapprox(nugget(γ), 0.05, atol=1e-3)

  # fix maximum parameters
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxrange=5.0u"m")
  @test isapprox(range(γ), 5.0u"m", atol=1e-3u"m")
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxsill=0.04)
  @test isapprox(sill(γ), 0.04, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxnugget=0.004)
  @test isapprox(nugget(γ), 0.004, atol=1e-3)

  # best fit is a Gaussian variogram
  γ = GeoStatsFunctions.fit(Variogram, g)
  @test γ isa GaussianVariogram
  @test isapprox(sill(γ), 0.054, atol=1e-3)
  γ = GeoStatsFunctions.fit([SphericalVariogram, GaussianVariogram], g)
  @test γ isa GaussianVariogram
  @test isapprox(sill(γ), 0.054, atol=1e-3)

  # make sure convenient methods work
  γ₁ = GeoStatsFunctions.fit(GaussianVariogram, g, h -> 1 / h)
  γ₂ = GeoStatsFunctions.fit(Variogram, g, h -> 1 / h)
  @test sill(γ₁) > 0
  @test sill(γ₂) > 0

  # unitful types
  img = readdlm(joinpath(datadir, "WalkerLake.txt"))
  d = georef((; Z=img * u"K"))
  g = EmpiricalVariogram(d, :Z, maxlag=15.0)
  γ = GeoStatsFunctions.fit(Variogram, g)
  @test unit(sill(γ)) == u"K^2"
  @test unit(nugget(γ)) == u"K^2"
  γ = GeoStatsFunctions.fit(Variogram, g, range=12.0u"m")
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")
  γ = GeoStatsFunctions.fit(Variogram, g, sill=0.06u"K^2")
  @test isapprox(sill(γ), 0.06u"K^2", atol=1e-3u"K^2")
  γ = GeoStatsFunctions.fit(Variogram, g, nugget=0.02u"K^2")
  @test isapprox(nugget(γ), 0.02u"K^2", atol=1e-3u"K^2")
  γ = GeoStatsFunctions.fit(Variogram, g, maxrange=6.0u"m")
  @test isapprox(range(γ), 6.0u"m", atol=1e-3u"m")
  γ = GeoStatsFunctions.fit(Variogram, g, sill=0.03u"K^2")
  @test isapprox(sill(γ), 0.03u"K^2", atol=1e-3u"K^2")
  γ = GeoStatsFunctions.fit(Variogram, g, nugget=0.01u"K^2")
  @test isapprox(nugget(γ), 0.01u"K^2", atol=1e-3u"K^2")

  # power variograms
  sₜ = 1.00
  nₜ = 0.20
  eₜ = 1.50
  γₜ = PowerVariogram(sₜ, nₜ, eₜ)
  xs = collect(0.0:0.1:10.0)u"m"
  ys = γₜ.(xs)
  ns = rand(1000:5000, length(xs))
  g = EmpiricalVariogram(ns, xs, ys, Euclidean(), :matheron)
  γ = GeoStatsFunctions.fit(PowerVariogram, g)
  @test isapprox(γ.nugget, nₜ, atol=1e-3)
  @test isapprox(γ.scaling, sₜ, atol=1e-3)
  @test isapprox(γ.exponent, eₜ, atol=1e-3)

  # test different settings
  sₜ = 6.54
  nₜ = 1.45
  eₜ = 0.64
  γₜ = PowerVariogram(sₜ, nₜ, eₜ)
  xs = collect(0.0:10.0:200.0)u"m"
  ys = γₜ.(xs)
  ns = rand(100:1000, length(xs))
  g = EmpiricalVariogram(ns, xs, ys, Euclidean(), :matheron)
  γ = GeoStatsFunctions.fit(PowerVariogram, g)
  @test isapprox(γ.nugget, nₜ, atol=1e-3)
  @test isapprox(γ.scaling, sₜ, atol=1e-3)
  @test isapprox(γ.exponent, eₜ, atol=1e-3)
end

@testset "Fitting transiograms" begin
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  t = EmpiricalTransiogram(gtb, "FACIES", maxlag=20, nlags=20)
  τ = GeoStatsFunctions.fit(PiecewiseLinearTransiogram, t)
  @test τ(0.0u"m") == I(5)
  @test all(x -> 0 < x < 1, τ(5.0u"m"))
  @test all(allequal, eachcol(τ(100.0u"m")))
end
