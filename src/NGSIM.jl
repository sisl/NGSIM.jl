__precompile__(true)

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
    load_trajdata,
    get_corresponding_roadway,
    convert_raw_ngsim_to_trajdatas


include("roadway.jl")
include("trajdata.jl")

end # module
