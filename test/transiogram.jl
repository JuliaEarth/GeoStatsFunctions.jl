@testset "Transiogram" begin
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = ExponentialTransiogram(R)
  @test t isa ExponentialTransiogram
  @test ratematrix(t) isa StaticMatrix
  @test range(t) == maximum(1 ./ -diag(R))

  # invalid transition rate matrix
  A = rand(3, 2)
  R = A ./ sum(A, dims=2)
  t = @test_throws ArgumentError ExponentialTransiogram(R)
end
