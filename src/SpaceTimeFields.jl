module SpaceTimeFields

  using UnPack
  using PolyLog

  import Base: +, -, *, show

  export EldredgeRamp, EldredgeRampIntegral, ColoniusRamp, Sinusoid, Gaussian, DGaussian, ConstantProfile, d_dt
  export EmptySpatialField, SpatialGaussian
  export AbstractSpatialField, Abstract1DProfile, SpatialField, SpatialTemporalField
  export FunctionProfile


  include("1dprofiles.jl")
  include("spatialfields.jl")


end
