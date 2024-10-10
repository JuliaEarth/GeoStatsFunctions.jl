@testset "Matrices" begin
  # helper function for tests
  function testratematrix(R)
    m, n = size(R)
    for i in 1:m, j in 1:n
      if i == j
        @test R[i, j] < 0 / u"m"
      else
        @test 0 / u"m" ≤ R[i, j] ≤ -R[i, i]
      end
    end
  end

  # synthetic geotable
  dom = CylindricalTrajectory(rand(Point, 100))
  tab = (; z=rand(1:3, 100))
  gtb = georef(tab, dom)

  # transition count matrix
  C = GeoStatsFunctions.countmatrix(gtb, "z")
  @test size(C) == (3, 3)
  @test eltype(C) == Int

  # transition probability matrix
  T = GeoStatsFunctions.probmatrix(gtb, "z")
  @test size(T) == (3, 3)
  @test eltype(T) == Float64

  # transition rate matrix
  R = GeoStatsFunctions.ratematrix(gtb, "z")
  @test size(R) == (3, 3)
  @test eltype(R) == typeof(1 / 1u"m")

  # basic properties of transition rate matrix
  csv = CSV.File(joinpath(datadir, "facies5.csv"))
  gtb = georef(csv, ("X", "Y", "Z"))
  R = GeoStatsFunctions.ratematrix(gtb, "FACIES")
  @test size(R) == (5, 5)
  testratematrix(R)

  # basic properties of transition rate matrix
  csv = CSV.File(joinpath(datadir, "facies15.csv"))
  gtb = georef(csv, ("X", "Y"))
  R = GeoStatsFunctions.ratematrix(gtb, "FACIES")
  @test size(R) == (15, 15)
  testratematrix(R)

  # base transition rate matrix
  R = GeoStatsFunctions.baseratematrix([1.0, 2.0, 3.0]u"m", [0.2, 0.5, 0.3])
  @test R ==
        [
    -1/1.0 0.5 / (1 - 0.2)/1.0 0.3 / (1 - 0.2)/1.0
    0.2 / (1 - 0.5)/2.0 -1/2.0 0.3 / (1 - 0.5)/2.0
    0.2 / (1 - 0.3)/3.0 0.5 / (1 - 0.3)/3.0 -1/3.0
  ] * u"m^-1"
end
