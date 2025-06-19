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

  # number of variates
  @test all(nvariates.(ts) .== 2)

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
  @test nvariates(t) == 3
  @test meanlengths(t) == (1.0u"m", 2.0u"m", 3.0u"m")
  @test all(proportions(t) .≈ (0.12403100775193801, 0.38759689922480606, 0.48837209302325607))
  @test range(GeoStatsFunctions.scale(t, 2)) == 6.0u"m"

  # random transition rate matrix
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = MatrixExponentialTransiogram(R)
  @test t isa MatrixExponentialTransiogram
  @test nvariates(t) == 3
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
  t = MatrixExponentialTransiogram(
    ranges=(2.0, 1.0),
    lengths=(1.0u"m", 2.0u"m", 3.0u"m"),
    proportions=(1 / 3, 1 / 3, 1 / 3)
  )
  @test meanlengths(t) == (2.0u"m", 4.0u"m", 6.0u"m")
end

@testset "PiecewiseLinearTransiogram" begin
  # basic tests
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  t = EmpiricalTransiogram(gtb, "FACIES", maxlag=20, nlags=20)
  τ = GeoStatsFunctions.fit(PiecewiseLinearTransiogram, t)
  @test range(GeoStatsFunctions.scale(τ, 2)) == 2u"m"
  @test nvariates(τ) == 5
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
  Y = [[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0], [0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5], [0.7 0.2 0.1; 0.7 0.2 0.1; 0.7 0.2 0.1]]
  τ = PiecewiseLinearTransiogram(h, Y)
  @test all(proportions(τ) .≈ (0.7, 0.2, 0.1))
  @test τ(100) ≈ [0.7 0.2 0.1; 0.7 0.2 0.1; 0.7 0.2 0.1]
end
