module SpaceTimeFields

  using UnPack

  import Base: +, -, *, show

  export EldredgeRamp, ColoniusRamp, Sinusoid, Gaussian, DGaussian, ConstantProfile, d_dt
  export EmptySpatialField, SpatialGaussian
  export AbstractSpatialField, Abstract1DProfile, SpatialField, SpatialTemporalField


  include("1dprofiles.jl")
  include("spatialfields.jl")


end
