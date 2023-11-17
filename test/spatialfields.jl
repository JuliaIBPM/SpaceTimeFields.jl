using Random

@testset "Profiles" begin

σ = rand()
g = Gaussian(σ,1)
@test g(σ) ≈ exp(-1)/sqrt(π*σ^2)

@test SpaceTimeFields.radius(g) == σ
@test SpaceTimeFields.strength(g) == 1

σ = 0.2
A = 1
dg = DGaussian(σ,A)
x = 0.1
@test dg(x) == -2*A*x/sqrt(π)/σ^3*exp(-x^2/σ^2)

a = 11
g = EldredgeRamp(a)
f = EldredgeRampIntegral(a)
g2 = d_dt(f)
x = rand()
@test g2(x) ≈ g(x)

end


@testset "SpatialFields" begin

σ = 0.2
A = 1

g = EmptySpatialField()
@test g(rand(2)...) == 0.0

g = SpatialGaussian(σ,σ,0.0,0.5,1)
@test g(0,0.5+σ) ≈ g(σ,0.5) ≈ g(-σ,0.5) ≈ g(0,0.5-σ) ≈ A*exp(-1/2)/(2π*σ^2)

u = 1
v = 0
gc = SpatialGaussian(σ,σ,0.0,0.5,1,u,v)
t = 1
@test gc(u*t,0.5+σ+v*t,t) ≈ gc(σ+u*t,0.5+v*t,t) ≈ gc(-σ+u*t,0.5+v*t,t) ≈ gc(0+u*t,0.5-σ+v*t,t) ≈ A*exp(-1/2)/(2π*σ^2)



g = EmptySpatialField()
for x in [-0.5,0,0.5], y in [-0.5,0,0.5]
  g += SpatialGaussian(0.2,0.5,x,y,1)
end
@test g(0.5,0) ≈ 3.6769641226297485

g2 = -g
@test g2(0.5,0) ≈ -3.6769641226297485



end

@testset "Derivatives of gaussians" begin


  σx = 0.5
  σy = 0.5
  x0 = 0
  y0 = 0
  A = 1
  dgaussx = SpatialGaussian(σx,σy,x0,y0,A,derivdir=1)
  dgaussy = SpatialGaussian(σx,σy,x0,y0,A,derivdir=2)
  x = 0.1
  y = 0.2
  @test dgaussx(0.1,0.2) == dgaussy(0.2,0.1) ≈ -A*x/(2π*σx^3*σy)*exp(-0.5*x^2/σx^2)*exp(-0.5*y^2/σy^2)

end
