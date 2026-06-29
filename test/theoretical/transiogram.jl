@testset "Transiogram" begin
  rng = StableRNG(123)
  h = range(0, stop=10, length=50)
  x, y = rand(rng, Point), rand(rng, Point)

  # simple transiogram models
  ts = [LinearTransiogram(), GaussianTransiogram(), SphericalTransiogram(), ExponentialTransiogram()]

  # check stationarity
  @test all(isstationary, ts)

  # check anisotropy
  @test all(isisotropic, ts)

  # check symmetry
  @test all(!issymmetric, ts)

  # check bandness
  @test all(isbanded, ts)

  # scale metric ball
  for t in ts
    t′ = GeoStatsFunctions.scale(t, 2)
    @test range(t′) == 2.0u"m"
  end

  # number of variables
  @test all(nvariables.(ts) .== 2)

  # mean lengths
  for t in ts
    @test meanlengths(t) == (1.0u"m", 1.0u"m")
  end

  # proportions
  for t in ts
    @test proportions(t) == (0.5, 0.5)
  end

  # identity matrix at lag zero
  for t in ts
    @test isapprox(t(0.0), [1.0 0.0; 0.0 1.0], atol=1e-5)
  end

  # converge to proportions
  for t in ts
    @test isapprox(t(1000.0), [0.5 0.5; 0.5 0.5], atol=1e-5)
  end

  # non-uniform proportions
  t = LinearTransiogram(proportions=(0.7, 0.2, 0.1))
  @test proportions(t) == (0.7, 0.2, 0.1)
  @test isapprox(t(1000.0), repeat([0.7 0.2 0.1], 3, 1), atol=1e-5)
  t = GaussianTransiogram(proportions=(0.7, 0.2, 0.1))
  @test proportions(t) == (0.7, 0.2, 0.1)
  @test isapprox(t(1000.0), repeat([0.7 0.2 0.1], 3, 1), atol=1e-5)

  # pairwise evaluation
  ps = [Point(0, 0), Point(1, 1), Point(2, 2)]
  for t in ts
    T = GeoStatsFunctions.pairwise(t, ps)
    @test size(T) == (6, 6)
  end
end

@testset "MatrixExponentialTransiogram" begin
  # base transition rate matrix
  R = GeoStatsFunctions.baseratematrix((1.0u"m", 2.0u"m", 3.0u"m"), (0.2, 0.5, 0.3))
  @test R ==
        [
    -1/1.0 0.5 / (1 - 0.2)/1.0 0.3 / (1 - 0.2)/1.0
    0.2 / (1 - 0.5)/2.0 -1/2.0 0.3 / (1 - 0.5)/2.0
    0.2 / (1 - 0.3)/3.0 0.5 / (1 - 0.3)/3.0 -1/3.0
  ] * u"m^-1"

  # corresponding exponential transiogram
  t = MatrixExponentialTransiogram(lengths=(1.0u"m", 2.0u"m", 3.0u"m"), proportions=(0.2, 0.5, 0.3))
  @test t isa MatrixExponentialTransiogram
  @test range(t) == 3.0u"m"
  @test nvariables(t) == 3
  @test meanlengths(t) == (1.0u"m", 2.0u"m", 3.0u"m")
  @test all(proportions(t) .≈ (0.12403100775193801, 0.38759689922480606, 0.48837209302325607))
  @test range(GeoStatsFunctions.scale(t, 2)) == 6.0u"m"

  # random transition rate matrix
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = MatrixExponentialTransiogram(R)
  @test t isa MatrixExponentialTransiogram
  @test nvariables(t) == 3
  @test meanlengths(t) == Tuple((1 ./ -diag(R)) * u"m")

  # invalid transition rate matrix
  A = rand(3, 2)
  R = A ./ sum(A, dims=2)
  t = @test_throws ArgumentError MatrixExponentialTransiogram(R)

  # pairwise evaluation
  ps = [Point(0, 0), Point(1, 1), Point(2, 2)]
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = MatrixExponentialTransiogram(R)
  T = GeoStatsFunctions.pairwise(t, ps)
  @test size(T) == (9, 9)

  # mean lengths with non-unit range
  t = MatrixExponentialTransiogram(range=3.0, lengths=(1.0u"m", 2.0u"m", 3.0u"m"), proportions=(1 / 3, 1 / 3, 1 / 3))
  @test meanlengths(t) == (3.0u"m", 6.0u"m", 9.0u"m")
  @test range(t) == 9.0u"m"
  t = MatrixExponentialTransiogram(
    ranges=(2.0, 1.0),
    lengths=(1.0u"m", 2.0u"m", 3.0u"m"),
    proportions=(1 / 3, 1 / 3, 1 / 3)
  )
  @test meanlengths(t) == (2.0u"m", 4.0u"m", 6.0u"m")
  @test range(t) == 6.0u"m"
end

@testset "PiecewiseLinearTransiogram" begin
  # basic tests
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  t = EmpiricalTransiogram(gtb, "FACIES", maxlag=20, nlags=20)
  τ = GeoStatsFunctions.fit(PiecewiseLinearTransiogram, t)
  @test range(GeoStatsFunctions.scale(τ, 2)) == 2range(τ)
  @test nvariables(τ) == 5
  @test τ(0.0u"m") == I(5)
  @test all(x -> 0 < x < 1, τ(5.0u"m"))
  @test all(allequal, eachcol(τ(100.0u"m")))

  # pairwise evaluation
  ps = [Point(0, 0, 0), Point(1, 1, 1), Point(2, 2, 2)]
  T = GeoStatsFunctions.pairwise(τ, ps)
  @test size(T) == (15, 15)

  # uniform proportions from zero off-diagonals
  h = [0.5, 1.5, 2.5]
  Y = [[1.0 0.0; 0.0 1.0], [0.5 0.5; 0.5 0.5], [0.0 1.0; 1.0 0.0]]
  τ = PiecewiseLinearTransiogram(h, Y)
  @test proportions(τ) == (0.5, 0.5)
  @test τ(100) == [0.5 0.5; 0.5 0.5]

  # non-uniform proportions
  h = [0.5, 1.5, 2.5]
  Y = [
    [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0],
    [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5],
    [0.7 0.2 0.1; 0.7 0.2 0.1; 0.7 0.2 0.1]
  ]
  τ = PiecewiseLinearTransiogram(h, Y)
  @test all(proportions(τ) .≈ (0.7, 0.2, 0.1))
  @test τ(100) ≈ [0.7 0.2 0.1; 0.7 0.2 0.1; 0.7 0.2 0.1]

  # custom ranges
  h = [0.5, 1.5, 2.5]
  Y = [[1.0 0.0; 0.0 1.0], [0.5 0.5; 0.5 0.5], [0.0 1.0; 1.0 0.0]]
  τ = PiecewiseLinearTransiogram(h, Y, ranges=(3.0, 2.0, 1.0))
  @test !isisotropic(τ)
  @test meanlengths(τ) == (7.5u"m", 7.5u"m")
  @test range(τ) == 7.5u"m"
  @test τ(Point(0, 0, 0), Point(1, 0, 0)) ≈ [1.0 0.0; 0.0 1.0]
  @test τ(Point(0, 0, 0), Point(0, 1, 0)) ≈ [1.0 0.0; 0.0 1.0]
  @test τ(Point(0, 0, 0), Point(0, 0, 1)) ≈ [0.75 0.25; 0.25 0.75]
end

@testset "CarleTransiogram" begin
  # anisotropic model from lengths and proportions
  l = (3.0u"m", 2.0u"m", 1.0u"m")
  p = (1 / 3, 1 / 3, 1 / 3)
  Rx = GeoStatsFunctions.baseratematrix(10 .* l, p)
  Ry = GeoStatsFunctions.baseratematrix(10 .* l, p)
  Rz = GeoStatsFunctions.baseratematrix(l, p)
  τ = CarleTransiogram(Rx, Ry, Rz)
  @test !isisotropic(τ)
  @test metricball(τ) == MetricBall(1.0u"m")
  @test meanlengths(τ) == 10 .* l
  @test range(τ) == 30.0u"m"
  @test τ(Point(0, 0, 0), Point(0, 0, 0)) == 1.0 * I(3)
  @test τ(Point(0, 0, 0), Point(1, 0, 0)) == τ(Point(0, 0, 0), Point(0, 1, 0))
  @test τ(Point(0, 0, 0), Point(1, 0, 0)) != τ(Point(0, 0, 0), Point(0, 0, 1))

  # isotropic models from single rate matrix
  τ = CarleTransiogram(Rx, Rx, Rx)
  @test isisotropic(τ)
  @test metricball(τ) == MetricBall(1.0u"m")
  @test meanlengths(τ) == 10 .* l
  @test range(τ) == 30.0u"m"
  τ = CarleTransiogram(Rz, Rz, Rz)
  @test isisotropic(τ)
  @test metricball(τ) == MetricBall(1.0u"m")
  @test meanlengths(τ) == l
  @test range(τ) == 3.0u"m"

  # anisotropic 2D model
  τ = CarleTransiogram(Rx, Rz)
  @test !isisotropic(τ)
  @test metricball(τ) == MetricBall(1.0u"m")
  @test meanlengths(τ) == 10 .* l
  @test range(τ) == 30.0u"m"
  @test τ(Point(0, 0), Point(0, 0)) == 1.0 * I(3)
  @test τ(Point(0, 0), Point(1, 0)) != τ(Point(0, 0), Point(0, 1))

  # default 3D isotropic model
  τ = CarleTransiogram()
  @test isisotropic(τ)
  @test meanlengths(τ) == (1.0u"m", 1.0u"m")
  @test proportions(τ) == (0.5, 0.5)

  # proportions must approximately sum up to one
  Rx = [
    -0.00477783 6.16608e-6 3.02809e-6 0.00185648 0.000770765 5.72372e-5 0.00208416
    2.34406e-5 -0.00927531 5.87122e-6 0.00359956 0.00149445 0.000110978 0.00404101
    8.75265e-6 4.46415e-6 -0.00346564 0.00134406 0.000558023 4.14389e-5 0.0015089
    2.66369e-7 1.35857e-7 6.6718e-8 -6.46327e-5 1.69823e-5 1.26111e-6 4.59203e-5
    3.96257e-7 2.02105e-7 9.92512e-8 6.08495e-5 -0.000131735 1.87606e-6 6.83121e-5
    1.08445e-6 5.53108e-7 2.71625e-7 0.000166529 6.9139e-5 -0.00042453 0.000186953
    4.89058e-8 2.49436e-8 1.22495e-8 7.51e-6 3.11797e-6 2.31542e-7 -1.09456e-5
  ]
  Ry = Rx
  Rz = [
    -0.0477783 6.16608e-5 3.02809e-5 0.0185648 0.00770765 0.000572372 0.0208416
    0.000234406 -0.0927531 5.87122e-5 0.0359956 0.0149445 0.00110978 0.0404101
    8.75265e-5 4.46415e-5 -0.0346564 0.0134406 0.00558023 0.000414389 0.015089
    2.66369e-6 1.35857e-6 6.6718e-7 -0.000646327 0.000169823 1.26111e-5 0.000459203
    3.96257e-6 2.02105e-6 9.92512e-7 0.000608495 -0.00131735 1.87606e-5 0.000683121
    1.08445e-5 5.53108e-6 2.71625e-6 0.00166529 0.00069139 -0.0042453 0.00186953
    4.89058e-7 2.49436e-7 1.22495e-7 7.51e-5 3.11797e-5 2.31542e-6 -0.000109456
  ]
  τ = CarleTransiogram(Rx, Ry, Rz)
  @test sum(proportions(τ)) ≈ 1 atol = 1e-2

  # meanlenghts with single rate matrix
  l = (3.0u"m", 2.0u"m", 1.0u"m")
  p = (1 / 3, 1 / 3, 1 / 3)
  R = GeoStatsFunctions.baseratematrix(l, p)
  τ = CarleTransiogram(R)
  @test meanlengths(τ) == (3.0u"m", 2.0u"m", 1.0u"m")
end
