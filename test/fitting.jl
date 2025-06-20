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

  # fix range without units
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, range=12.0)
  @test isapprox(range(γ), 12.0u"m", atol=1e-3u"m")

  # fix maximum parameters
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxrange=5.0u"m")
  @test isapprox(range(γ), 5.0u"m", atol=1e-3u"m")
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxsill=0.04)
  @test isapprox(sill(γ), 0.04, atol=1e-3)
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxnugget=0.004)
  @test isapprox(nugget(γ), 0.004, atol=1e-3)

  # fix maximum range without units
  γ = GeoStatsFunctions.fit(GaussianVariogram, g, maxrange=5.0)
  @test isapprox(range(γ), 5.0u"m", atol=1e-3u"m")

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
  γₜ = PowerVariogram(; scaling=sₜ, nugget=nₜ, exponent=eₜ)
  xs = collect(0.0:0.1:10.0)u"m"
  ys = γₜ.(xs)
  ns = rand(1000:5000, length(xs))
  g = EmpiricalVariogram(ns, xs, ys, Euclidean(), :matheron)
  γ = GeoStatsFunctions.fit(PowerVariogram, g)
  @test isapprox(γ.nugget, nₜ, atol=1e-3)
  @test isapprox(γ.scaling, sₜ, atol=1e-3)
  @test isapprox(γ.exponent, eₜ, atol=1e-3)
  sₜ = 6.54
  nₜ = 1.45
  eₜ = 0.64
  γₜ = PowerVariogram(; scaling=sₜ, nugget=nₜ, exponent=eₜ)
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

  # all fits lead to similar proportions
  τ₁ = GeoStatsFunctions.fit(LinearTransiogram, t)
  τ₂ = GeoStatsFunctions.fit(SphericalTransiogram, t)
  τ₃ = GeoStatsFunctions.fit(GaussianTransiogram, t)
  τ₄ = GeoStatsFunctions.fit(ExponentialTransiogram, t)
  ps = (0.30, 0.12, 0.12, 0.20, 0.25)
  @test all(isapprox.(proportions(τ₁), ps, atol=1e-2))
  @test all(isapprox.(proportions(τ₂), ps, atol=1e-2))
  @test all(isapprox.(proportions(τ₃), ps, atol=1e-2))
  @test all(isapprox.(proportions(τ₄), ps, atol=1e-2))

  # fix parameters
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, range=12.0u"m")
  @test isapprox(range(τ), 12.0u"m", atol=1e-3u"m")
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, proportions=ntuple(i -> 1 / 5, 5))
  @test all(isapprox.(proportions(τ), 1 / 5, atol=1e-3))

  # fix range without units
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, range=12.0)
  @test isapprox(range(τ), 12.0u"m", atol=1e-3u"m")

  # fix maximum parameters
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, maxrange=5.0u"m")
  @test isapprox(range(τ), 5.0u"m", atol=1e-3u"m")
  ps = (0.2, 0.1, 0.1, 0.3, 0.25)
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, maxproportions=ps)
  @test all(isapprox.(proportions(τ), ps, atol=1e-1))

  # fix maximum range without units
  τ = GeoStatsFunctions.fit(SphericalTransiogram, t, maxrange=5.0)
  @test isapprox(range(τ), 5.0u"m", atol=1e-3u"m")

  # matrix-exponential
  τ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t)
  ps = (0.30, 0.12, 0.12, 0.20, 0.25)
  @test all(isapprox.(proportions(τ), ps, atol=3e-2))
  τ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t, range=0.8u"m")
  @test isapprox(radius(metricball((τ))), 0.8u"m", atol=1e-3u"m")
  τ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t, range=0.8)
  @test isapprox(radius(metricball((τ))), 0.8u"m", atol=1e-3u"m")
  ls = (7.0u"m", 3.0u"m", 2.0u"m", 7.0u"m", 14.0u"m")
  τ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t, range=1.0u"m", maxlengths=ls)
  @test all(isapprox.(meanlengths(τ)[1:2], ls[1:2], atol=1e-3u"m"))
  @test all(meanlengths(τ)[3:5] .≤ ls[3:5])
  ps = (0.32, 0.1, 0.1, 0.3, 0.31)
  τ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t, maxproportions=ps)
  @test all(proportions(τ) .≤ ps)

  # piecewise linear
  τ = GeoStatsFunctions.fit(PiecewiseLinearTransiogram, t)
  @test τ(0.0u"m") == I(5)
  @test all(x -> 0 < x < 1, τ(5.0u"m"))
  @test all(allequal, eachcol(τ(100.0u"m")))

  # best fit is a matrix-exponential transiogram
  ps = (0.30, 0.12, 0.12, 0.20, 0.25)
  τ = GeoStatsFunctions.fit(Transiogram, t)
  @test τ isa MatrixExponentialTransiogram
  @test all(isapprox.(proportions(τ), ps, atol=3e-2))
  τ = GeoStatsFunctions.fit([SphericalTransiogram, MatrixExponentialTransiogram], t)
  @test τ isa MatrixExponentialTransiogram
  @test all(isapprox.(proportions(τ), ps, atol=3e-2))

  # make sure convenient methods work
  ps = (0.30, 0.12, 0.12, 0.20, 0.25)
  τ₁ = GeoStatsFunctions.fit(MatrixExponentialTransiogram, t, h -> 1 / h)
  τ₂ = GeoStatsFunctions.fit(Transiogram, t, h -> 1 / h)
  @test all(isapprox.(proportions(τ₁), ps, atol=3e-2))
  @test all(isapprox.(proportions(τ₂), ps, atol=3e-2))
end
