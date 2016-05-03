module NGSIM

using Vec
using DataFrames
using Distributions

export
    CurvePt,
    Roadway,
    ROADWAY_80,
    ROADWAY_101,

    TrajdataRaw,
    Trajdata,
    Frenet,
    VehicleSystem,
    VehicleState,
    Vehicle,
    FilterTrajectoryResult,


    read_roadway,
    get_roadway_for_trajdata,

    project_to_lane,
    project_to_closest_lane,
    curve_at,
    move_extind_along,
    get_neighbor_laneid_left,
    get_neighbor_laneid_right,

    SMOOTHING_WIDTH_POS,
    CLASS_MOTORCYCLE,
    CLASS_AUTOMOBILE,
    CLASS_TRUCKORBUS,

    get_vehicle,
    get_vehicle!,
    get_vehicles,
    get_vehiclestate,

    get_footpoint,
    get_center,
    get_frame_range,
    get_index_of_first_vehicle_with_id,
    nframes,
    carsinframe,
    carid_set,
    nth_carid,
    first_carid,
    iscarinframe,
    car_df_index,

    get_state_list,
    get_state_list_global,
    get_state_list_frenet,

    get_neighbor_index_fore,
    get_neighbor_index_rear,
    get_headway_dist_between,
    get_headway_time_between,

    filter_trajectory!,
    symmetric_exponential_moving_average!,
    load_trajdata_raw,
    input_path_to_extracted_trajdata_csv

include("roadway.jl")
include("trajdata.jl")

end # module
