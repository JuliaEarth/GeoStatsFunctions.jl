@testset "EmpiricalTransiogram" begin
  # diagonal ordinates are non-increasing
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  t = EmpiricalTransiogram(gtb, "FACIES", nlags=20)
  for i in 1:5
    ys = t.ordinates[i, i]
    @test ys[1] > ys[10] > ys[20]
  end

  # print methods
  rng = StableRNG(123)
  d = georef((; z=rand(rng, 1:10, 1000)), rand(rng, Point, 1000))
  t = EmpiricalTransiogram(d, :z)
  @test sprint(show, t) == "EmpiricalTransiogram(distance: Euclidean(0.0), estimator: CarleEstimator(), npairs: 176100)"
  @test sprint(show, MIME"text/plain"(), t) == """
  EmpiricalTransiogram
  ├─ abscissas: [0.00249657 m, 0.00828093 m, 0.0128198 m, ..., 0.0875733 m, 0.0925317 m, 0.0974132 m]
  ├─ ordinates: 
  │  ├─ [0.0, 0.0, 0.0, ..., 0.0, 0.0, 0.107143]
  │  ├─ [0.0, 0.0, 0.0, ..., 0.0, 0.0869565, 0.103448]
  │  ├─ [0.0, 0.0, 0.0, ..., 0.130435, 0.0909091, 0.142857]
  │  ⋮
  │  ├─ [0.0, 0.0, 0.0, ..., 0.0555556, 0.166667, 0.0833333]
  │  ├─ [0.0, 0.0, 0.0, ..., 0.125, 0.208333, 0.1]
  │  └─ [0.0, 0.0, 0.0, ..., 0.130435, 0.0833333, 0.0909091]
  ├─ distance: Euclidean(0.0)
  ├─ estimator: CarleEstimator()
  └─ npairs: 176100"""
end
