using SpaceTimeFields
using Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "SpatialFields"
  include("spatialfields.jl")
end
