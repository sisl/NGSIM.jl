__precompile__(true)

module NGSIM

using AutomotiveDrivingModels
using DataFrames
using Distributions

export
    NGSIMRoadway,
    RoadwayInputParams,

    ROADWAY_80,
    ROADWAY_101,

    NGSIMTrajdata,
    # VehicleSystem,
    # FilterTrajectoryResult,

    NGSIM_TIMESTEP,
    NGSIM_TRAJDATA_PATHS,
    # TRAJDATA_PATHS,

    carsinframe,
    load_ngsim_trajdata,
    get_corresponding_roadway
    # filter_trajectory!,
    # symmetric_exponential_moving_average!,
    # load_trajdata,
    # convert_raw_ngsim_to_trajdatas


include("roadway.jl")
include("ngsim_trajdata.jl")
# include("trajdata.jl")

end # module
