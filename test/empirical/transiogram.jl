@testset "EmpiricalTransiogram" begin
  # diagonal ordinates are non-increasing
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  t = EmpiricalTransiogram(gtb, "FACIES", maxlag=20, nlags=20)
  for i in 1:5
    ys = t.ordinates[i, i]
    @test ys[1] > ys[10] > ys[20]
  end

  # same for directional transiograms
  t = DirectionalTransiogram((0.0, 0.0, 1.0), gtb, "FACIES", maxlag=20, nlags=20)
  for i in 1:5
    ys = t.ordinates[i, i]
    @test ys[1] > ys[10] > ys[20]
  end
end
