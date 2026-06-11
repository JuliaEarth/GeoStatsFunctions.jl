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

  # additional columns don't affect results
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  tab = (X=csv.X, Y=csv.Y, Z=csv.Z, FACIES=csv.FACIES, NAME="FACIES" .* string.(csv.FACIES))
  gtb1 = georef(csv, ("X", "Y", "Z"))
  gtb2 = georef(tab, ("X", "Y", "Z"))
  t1 = EmpiricalTransiogram(gtb1, "FACIES", maxlag=20, nlags=20)
  t2 = EmpiricalTransiogram(gtb2, "FACIES", maxlag=20, nlags=20)
  @test t1.ordinates == t2.ordinates
end
