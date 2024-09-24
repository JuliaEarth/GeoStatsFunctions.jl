@testset "Transiogram" begin
  A = rand(3, 3)
  R = A ./ sum(A, dims=2)
  t = ExponentialTransiogram(R)
  @test t isa ExponentialTransiogram
end
