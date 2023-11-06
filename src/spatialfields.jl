### Routines for constructing smooth functions that generate field data ###

using LinearAlgebra
using StaticArrays
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

## Basic spatial field ##
"""
    SpatialField(f::Function)

Creates a lazy instance of a spatial field from a function `f`, which
must have the signature `f(x,y)`.
"""
struct SpatialField{FT<:Function} <: AbstractSpatialField
  f :: FT
end
(field::SpatialField)(x,y) = field.f(x,y)
(field::SpatialField)(x,y,t) = field.f(x,y) # ignore the time argument

"""
    SpatialTemporalField(f::Function)

Creates a lazy instance of a spatial-temporal field from a function `f`, which
must have the signature `f(x,y,t)`.
"""
struct SpatialTemporalField{FT<:Function} <: AbstractSpatialField
  f :: FT
end
(field::SpatialTemporalField)(x,y,t) = field.f(x,y,t)


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
""" SpatialGaussian(::Real,::Real,::Real,::Real,::Real)

"""
    SpatialGaussian(Σ::Matrix,x0::Vector,A::Real[,derivdir=0])

Set up a spatial field in the form of a Gaussian centered at `x0[1],x0[2]` with
covariance matrix `Σ` and amplitude `A`. If the
optional parameter `deriv` is set to 1 or 2, then it returns the first
derivative of a Gaussian in that direction (`x` or `y`, respectively).
""" SpatialGaussian(::AbstractMatrix,::Vector,::Real)

"""
    SpatialGaussian(σ::Vector,x0::Real,y0::Real,α::Real,A::Real[,derivdir=0])

Set up a spatial field in the form of a Gaussian centered at `x0,y0` with
variances `σ[1],σ[2]` along the orthoogonal eigendirections, which are rotated by `α`
with respect to the coordinate system; and amplitude `A`. If the
optional parameter `deriv` is set to 1 or 2, then it returns the first
derivative of a Gaussian in that direction (`x` or `y`, respectively).
""" SpatialGaussian(::Vector,::Real,::Real,::Real,::Real)


struct SpatialGaussian{CT,D} <: AbstractSpatialField
  Σ :: SMatrix{2,2,Float64,4}
  Σinv :: SMatrix{2,2,Float64,4}
  x0 :: SVector{2,Float64}
  A :: Float64
  fact :: Float64
  u :: SVector{2,Float64}
  SpatialGaussian(Σ::SMatrix,x0::SVector,A::Real;derivdir=0) = new{false,derivdir}(Σ,inv(Σ),x0,A,A/(2π*sqrt(det(Σ))),SVector{2}(0.0,0.0))
  SpatialGaussian(Σ::SMatrix,x0::SVector,A::Real,u::SVector;derivdir=0) = new{true,derivdir}(Σ,inv(Σ),x0,A,A/(2π*sqrt(det(Σ))),u)
  #SpatialGaussian(gx::Abstract1DProfile,gy::Abstract1DProfile,A,u,v) = new{true,typeof(gx),typeof(gy)}(gx,gy,A,u,v)
  #SpatialGaussian(gx::Abstract1DProfile,gy::Abstract1DProfile,A) = new{false,typeof(gx),typeof(gy)}(gx,gy,A,0.0,0.0)
end


function SpatialGaussian(σ2::Vector,x0::Real,y0::Real,α::Real,A::Real;kwargs...)
    R = _rotation_matrix(α)
    Σ = R*SMatrix{2,2}(Diagonal(σ2))*R'
    x0v = SVector{2}(x0,y0)
    SpatialGaussian(Σ,x0v,A;kwargs...)
end

function SpatialGaussian(σ2::Vector,x0::Real,y0::Real,α::Real,A::Real,u::Real,v::Real;kwargs...)
    R = _rotation_matrix(α)
    Σ = R*SMatrix{2,2}(Diagonal(σ2))*R'
    x0v = SVector{2}(x0,y0)
    uv = SVector{2}(u,v)
    SpatialGaussian(Σ,x0v,A,uv;kwargs...)
end

SpatialGaussian(Σ::AbstractMatrix,x0::Vector,a...;kwargs...) = SpatialGaussian(SMatrix{2,2}(Σ),SVector{2}(x0),a...;kwargs...)

SpatialGaussian(σx::Real,σy::Real,x0::Real,y0::Real,a...;kwargs...) = SpatialGaussian([σx^2,σy^2],x0,y0,0.0,a...;kwargs...)



function _rotation_matrix(θ::Real)
    cth, sth = cos(θ), sin(θ)
    R = @SMatrix[cth -sth; sth cth]
    return R
end

#=
SpatialGaussian(σx::Real,σy::Real,x0::Real,y0::Real,A::Real;deriv::Int=0) =
                _spatialdgaussian(σx,σy,x0,y0,A,Val(deriv))
SpatialGaussian(σx::Real,σy::Real,x0::Real,y0::Real,A::Real,u::Real,v::Real;deriv::Int=0) =
                _spatialdgaussian(σx,σy,x0,y0,A,u,v,Val(deriv))
=#

#=
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
=#

function (g::SpatialGaussian{C,0})(x,y) where {C}
    @unpack x0, Σinv, fact = g
    xv = SVector{2}(x,y)
    return _spatialgaussian(fact,xv-x0,Σinv)
end

function (g::SpatialGaussian{C,N})(x,y) where {C,N}
    @unpack x0, Σinv, fact = g
    xv = SVector{2}(x,y)
    dx = xv-x0
    dfact = -Σinv*dx
    return _spatialgaussian(fact,dx,Σinv)*dfact[N]
end

(g::SpatialGaussian{false,N})(x,y,t) where {N} = g(x,y)
(g::SpatialGaussian{true,N})(x,y,t) where {N} = g(x-g.u[1]*t,y-g.u[2]*t)

function (g::SpatialGaussian{C,N})(x,y) where {C,N}
    @unpack x0, Σinv, fact = g
    xv = SVector{2}(x,y)
    dx = xv-x0
    dfact = -Σinv*dx
    return _spatialgaussian(fact,dx,Σinv)*dfact[N]
end

_spatialgaussian(fact::Real,x::SVector,Σinv::SMatrix) = fact*exp(-0.5*x'*Σinv*x)

#=
(g::SpatialGaussian{GX,GY})(x,y) where {GX,GY} = g.gx(x)*g.gy(y)
# ignore the time argument if it is called with this...
(g::SpatialGaussian{false,GX,GY})(x,y,t) where {GX,GY} = g(x,y)
(g::SpatialGaussian{true,GX,GY})(x,y,t) where {GX,GY} = g.gx(x-g.u*t)*g.gy(y-g.v*t)
=#

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
