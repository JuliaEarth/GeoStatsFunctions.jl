@testset "MatrixExponentialTransiogram" begin
  # base transition rate matrix
  R = GeoStatsFunctions.baseratematrix([1.0, 2.0, 3.0]u"m", [0.2, 0.5, 0.3])
  @test R ==
        [
    -1/1.0 0.5 / (1 - 0.2)/1.0 0.3 / (1 - 0.2)/1.0
    0.2 / (1 - 0.5)/2.0 -1/2.0 0.3 / (1 - 0.5)/2.0
    0.2 / (1 - 0.3)/3.0 0.5 / (1 - 0.3)/3.0 -1/3.0
  ] * u"m^-1"

  # corresponding exponential transiogram
  t = MatrixExponentialTransiogram([1.0, 2.0, 3.0]u"m", [0.2, 0.5, 0.3])
  @test t isa MatrixExponentialTransiogram
  @test GeoStatsFunctions.ranges(t) == [1.0, 2.0, 3.0]u"m"
  @test range(t) == 3.0u"m"

  # random transition rate matrix
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = MatrixExponentialTransiogram(R)
  @test t isa MatrixExponentialTransiogram
  @test range(t) == maximum(1 ./ -diag(R))

  # invalid transition rate matrix
  A = rand(3, 2)
  R = A ./ sum(A, dims=2)
  t = @test_throws ArgumentError MatrixExponentialTransiogram(R)
end
