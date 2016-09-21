VERSION >= v"0.4.0-dev+6521" && __precompile__(true)


module NGSIM

using Compat
using AutomotiveDrivingModels
using DataFrames
using Distributions

export
    NGSIMRoadway,
    RoadwayInputParams,

    ROADWAY_80,
    ROADWAY_101,

    NGSIMTrajdata,
    VehicleSystem,
    FilterTrajectoryResult,

    TRAJDATA_PATHS,
    NGSIM_TIMESTEP,

    filter_trajectory!,
    symmetric_exponential_moving_average!,
    load_ngsim_trajdata,
    load_trajdata


include("roadway.jl")
include("trajdata.jl")

end # module
