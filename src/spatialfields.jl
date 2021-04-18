### Routines for constructing smooth functions that generate field data ###

abstract type AbstractSpatialField end

## Empty spatial field

"""
    EmptySpatialField()

Create a blank spatial field. This is primarily useful for initializing a
sum of spatial fields.

# Example
```jldoctest
julia> g = EmptySpatialField()
EmptySpatialField()

julia> g(2,3)
0.0
```
"""
struct EmptySpatialField <: AbstractSpatialField end
(g::EmptySpatialField)(a...) = Float64(0)


## Spatial Gaussian field ##

"""
    SpatialGaussian(σx,σy,x0,y0,A[,derivdir=0])

Set up a spatial field in the form of a Gaussian centered at `x0,y0` with
radii `σx` and `σy` in the respective directions and amplitude `A`. If the
optional parameter `deriv` is set to 1 or 2, then it returns the first
derivative of a Gaussian in that direction (`x` or `y`, respectively).

`SpatialGaussian(σx,σy,x0,y0,A,u,v[,derivdir=0])` generates a Gaussian
that convects at velocity `(u,v)`. It can be evaluated with an additional
argument for time.
"""
struct SpatialGaussian{CT,GX,GY} <: AbstractSpatialField
  gx :: GX
  gy :: GY
  A :: Float64
  u :: Float64
  v :: Float64
  SpatialGaussian(gx::Abstract1DProfile,gy::Abstract1DProfile,A,u,v) = new{true,typeof(gx),typeof(gy)}(gx,gy,A,u,v)
  SpatialGaussian(gx::Abstract1DProfile,gy::Abstract1DProfile,A) = new{false,typeof(gx),typeof(gy)}(gx,gy,A,0.0,0.0)
end

SpatialGaussian(σx::Real,σy::Real,x0::Real,y0::Real,A::Real;deriv::Int=0) =
                _spatialdgaussian(σx,σy,x0,y0,A,Val(deriv))
SpatialGaussian(σx::Real,σy::Real,x0::Real,y0::Real,A::Real,u::Real,v::Real;deriv::Int=0) =
                _spatialdgaussian(σx,σy,x0,y0,A,u,v,Val(deriv))


_spatialdgaussian(σx,σy,x0,y0,A,::Val{0}) = SpatialGaussian(Gaussian(σx,A) >> x0,
                                                            Gaussian(σy,1) >> y0,A)
_spatialdgaussian(σx,σy,x0,y0,A,::Val{1}) = SpatialGaussian(DGaussian(σx,A) >> x0,
                                                            Gaussian(σy,1) >> y0,A)
_spatialdgaussian(σx,σy,x0,y0,A,::Val{2}) = SpatialGaussian(Gaussian(σx,A) >> x0,
                                                            DGaussian(σy,1) >> y0,A)

_spatialdgaussian(σx,σy,x0,y0,A,u,v,::Val{0}) = SpatialGaussian(Gaussian(σx,A) >> x0,
                                                                Gaussian(σy,1) >> y0,A,u,v)
_spatialdgaussian(σx,σy,x0,y0,A,u,v,::Val{1}) = SpatialGaussian(DGaussian(σx,A) >> x0,
                                                                Gaussian(σy,1) >> y0,A,u,v)
_spatialdgaussian(σx,σy,x0,y0,A,u,v,::Val{2}) = SpatialGaussian(Gaussian(σx,A) >> x0,
                                                                DGaussian(σy,1) >> y0,A,u,v)


SpatialGaussian(σ,x0,y0,A;deriv::Int=0) = SpatialGaussian(σ,σ,x0,y0,A,deriv=deriv)


(g::SpatialGaussian{GX,GY})(x,y) where {CT,GX,GY} = g.gx(x)*g.gy(y)
# ignore the time argument if it is called with this...
(g::SpatialGaussian{false,GX,GY})(x,y,t) where {GX,GY} = g(x,y)
(g::SpatialGaussian{true,GX,GY})(x,y,t) where {GX,GY} = g.gx(x-g.u*t)*g.gy(y-g.v*t)


## Scaling spatial fields

struct ScaledField{N <: Real, P <: AbstractSpatialField} <: AbstractSpatialField
    s::N
    p::P
end
function show(io::IO, p::ScaledField)
    print(io, "$(p.s) × ($(p.p))")
end
s::Number * p::AbstractSpatialField = ScaledField(s, p)
-(p::AbstractSpatialField) = ScaledField(-1, p)

(p::ScaledField)(x,y) = p.s*p.p(x,y)
(p::ScaledField)(x,y,t) = p.s*p.p(x,y,t)


## Adding spatial fields together.

struct AddedFields{T <: Tuple} <: AbstractSpatialField
    ps::T
end
function show(io::IO, Σp::AddedFields)
    println(io, "AddedFields:")
    for p in Σp.ps
        println(io, "  $p")
    end
end

"""
    p₁::AbstractSpatialField + p₂::AbstractSpatialField

Add the fields so that `(p₁ + p₂)(x,y) = p₁(x,y) + p₂(x,y)`.

"""
+(p::AbstractSpatialField, Σp::AddedFields) = AddedFields((p, Σp.ps...))
+(Σp::AddedFields, p::AbstractSpatialField) = AddedFields((Σp.ps..., p))
function +(Σp₁::AddedFields, Σp₂::AddedFields)
    AddedFields((Σp₁..., Σp₂...))
end

-(p₁::AbstractSpatialField, p₂::AbstractSpatialField) = p₁ + (-p₂)

+(p::AbstractSpatialField...) = AddedFields(p)

# Evaluate at x,y , assuming that t = 0 for any time-varying member
function (Σp::AddedFields)(x,y)
    f = 0.0
    for p in Σp.ps
        f += p(x,y,0.0)
    end
    f
end

# Evaluate at x,y,t, ignoring t for any constant member
function (Σp::AddedFields)(x,y,t)
    f = 0.0
    for p in Σp.ps
        f += p(x,y,t)
    end
    f
end
