using SpaceTimeFields
using Test

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "SpatialFields"
  include("spatialfields.jl")
end

if GROUP == "All" || GROUP == "Aqua"
  include("aqua.jl")
end