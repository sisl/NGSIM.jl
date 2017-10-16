using Base.Test
using NGSIM
using NBInclude

let
    roadway_80 = open(io->read_dxf(io, Roadway, dist_threshold_lane_connect=2.0), Pkg.dir("NGSIM", "data", "ngsim_80.dxf"), "r")
end

# nbinclude(joinpath(dirname(@__FILE__), "..", "jnotebooks", "Demo.ipynb"))
