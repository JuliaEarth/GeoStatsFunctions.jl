@testset "Matrices" begin
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
end
