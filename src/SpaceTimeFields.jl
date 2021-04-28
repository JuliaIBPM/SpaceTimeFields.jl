module SpaceTimeFields

  import Base: +, -, *, show

  export EldredgeRamp, ColoniusRamp, Sinusoid, Gaussian, DGaussian, ConstantProfile, d_dt
  export EmptySpatialField, SpatialGaussian


  include("1dprofiles.jl")
  include("spatialfields.jl")


end
