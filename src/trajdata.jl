const NGSIM_TIMESTEP = 0.1 # [sec]
const SMOOTHING_WIDTH_POS = 0.5 # [s]

const CLASS_MOTORCYCLE = 1
const CLASS_AUTOMOBILE = 2
const CLASS_TRUCKORBUS = 3

include(Pkg.dir("NGSIM", "src", "trajectory_smoothing.jl"))

# from "Estimating Acceleration and Lane-Changing
#       Dynamics Based on NGSIM Trajectory Data"
function symmetric_exponential_moving_average(
    arr :: Vector{Float64},
    T   :: Float64; # smoothing width [s]
    dt  :: Float64 = 0.1 # sampling period [s]
    )

    Δ = T / dt

    N = length(arr)
    retval = Array(Float64, N)

    for i = 1 : N

        Z = 0.0
        x = 0.0

        D = min(round(Int, 3Δ), i-1)
        if i+D > N
            D = N-i
        end

        for k in (i-D):(i+D)
            e = exp(-abs(i-k)/Δ)
            Z += e
            x += arr[k] * e
        end

        retval[i] = x / Z
    end

    retval
end

const TRAJDATA_INPUT_PATHS = [
    "/media/tim/DATAPART1/Data/NGSIM/HW80/I-80-Main-Data/vehicle-trajectory-data/0400pm-0415pm/trajectories-0400-0415.txt",
    "/media/tim/DATAPART1/Data/NGSIM/HW80/I-80-Main-Data/vehicle-trajectory-data/0500pm-0515pm/trajectories-0500-0515.txt",
    "/media/tim/DATAPART1/Data/NGSIM/HW80/I-80-Main-Data/vehicle-trajectory-data/0515pm-0530pm/trajectories-0515-0530.txt",
    "/media/tim/DATAPART1/Data/NGSIM/HW101/US-101-Main-Data/vehicle-trajectory-data/0750am-0805am/trajectories-0750am-0805am.txt",
    "/media/tim/DATAPART1/Data/NGSIM/HW101/US-101-Main-Data/vehicle-trajectory-data/0805am-0820am/trajectories-0805am-0820am.txt",
    "/media/tim/DATAPART1/Data/NGSIM/HW101/US-101-Main-Data/vehicle-trajectory-data/0820am-0835am/trajectories-0820am-0835am.txt",
]

###############

immutable Frenet
    laneid::Int
    extind::Float64
    s::Float64 # distance along lane
    t::Float64 # lane offset, positive is to left
    ϕ::Float64 # lane relative heading
end
immutable VehicleState
    posG::VecSE2 # global
    posF::Frenet # (extind,t,ϕ)
    v::Float64

    VehicleState() = new(VecSE2(), Frenet(-1, NaN, NaN, NaN, NaN), NaN)
    VehicleState(posG::VecSE2, v::Float64) = new(posG, Frenet(-1, NaN, NaN, NaN, NaN), v)
    VehicleState(posG::VecSE2, posF::Frenet, v::Float64) = new(posG, posF, v)
end
type Vehicle
    id::Int
    class::Int # ∈ (1-motorcycle, 2-auto, 3-truck)
    length::Float64
    width::Float64
    state::VehicleState

    Vehicle() = new(0, 0, NaN, NaN, VehicleState())
end

get_footpoint(veh::Vehicle) = veh.state.posG + polar(veh.state.posF.t, veh.state.posG.θ-veh.state.posF.ϕ-π/2)
get_center(veh::Vehicle) = veh.state.posG + polar(veh.length/2, veh.state.posG.θ+π)

###############

type TrajdataRaw
    df         :: DataFrame
    car2start  :: Dict{Int, Int}         # maps carindex to starting index in the df
    frame2cars :: Dict{Int, Vector{Int}} # maps frame to list of carids in the scene
    roadway    :: Roadway

    function TrajdataRaw(input_path::AbstractString, roadway::Roadway)

        @assert(isfile(input_path))

        if splitext(input_path)[2] == ".txt" # txt is original
            df = readtable(input_path, separator=' ', header = false)
            col_names = [:id, :frame, :n_frames_in_dataset, :epoch, :local_x, :local_y, :global_x, :global_y, :length, :width, :class, :speed, :acc, :lane, :carind_front, :carind_rear, :dist_headway, :time_headway]
            for (i,name) in enumerate(col_names)
                rename!(df, symbol(@sprintf("x%d", i)), name)
            end

            df[:global_heading] = fill(NaN, nrow(df))
            df[:laneid] = fill(-1, nrow(df))
            df[:frenet_extind] = fill(NaN, nrow(df))
            df[:frenet_s] = fill(NaN, nrow(df))
            df[:frenet_t] = fill(NaN, nrow(df))
            df[:frenet_heading] = fill(NaN, nrow(df))
        else
            @assert(splitext(input_path)[2] == ".csv")
            df = readtable(input_path) # csv file is exported after extracting everything
        end

        car2start = Dict{Int, Int}()
        frame2cars = Dict{Int, Vector{Int}}()

        for (dfind, carid) in enumerate(df[:id])
            if !haskey(car2start, carid)
                car2start[carid] = dfind
            end

            frame = convert(Int, df[dfind, :frame])
            if !haskey(frame2cars, frame)
                frame2cars[frame] = [carid]
            else
                frame2cars[frame] = push!(frame2cars[frame], carid)
            end
        end

        new(df, car2start, frame2cars, roadway)
    end
end

nframes(trajdata::TrajdataRaw) = maximum(keys(trajdata.frame2cars))
carsinframe(trajdata::TrajdataRaw, frame::Int) = get(trajdata.frame2cars, frame, Int[]) # NOTE: memory allocation!
carid_set(trajdata::TrajdataRaw) = Set(keys(trajdata.car2start)) # NOTE: memory allocation!
nth_carid(trajdata::TrajdataRaw, frame::Int, n::Int) = trajdata.frame2cars[frame][n]
first_carid(trajdata::TrajdataRaw, frame::Int) = nth_carid(trajdata, frame, 1)
iscarinframe(trajdata::TrajdataRaw, carid::Int, frame::Int) = in(carid, carsinframe(trajdata, frame))

function get_vehicle!(
    veh::Vehicle,
    trajdata::TrajdataRaw,
    carid::Int,
    frame::Int,
    )

    dfind = car_df_index(trajdata, carid, frame)
    df = trajdata.df

    veh.id = carid
    veh.class  = df[dfind, :class ]
    veh.length = df[dfind, :length]
    veh.width  = df[dfind, :width ]

    posG = VecSE2(df[dfind, :global_x], df[dfind, :global_y], df[dfind, :global_heading])
    posF = Frenet(df[dfind, :laneid], df[dfind, :frenet_extind], df[dfind, :frenet_s], df[dfind, :frenet_t], df[dfind, :frenet_heading])
    speed = df[dfind, :speed]
    veh.state = VehicleState(posG, posF, speed)

    veh
end
function car_df_index(trajdata::TrajdataRaw, carid::Int, frame::Int)
    #=
    given frame and carid, find index of car in trajdata
    Returns 0 if it does not exist
    =#

    df = trajdata.df

    lo = trajdata.car2start[carid]
    framestart = df[lo, :frame]

    retval = 0

    if framestart == frame
        retval = lo
    elseif frame ≥ framestart
        retval = frame - framestart + lo
        n_frames = df[lo, :n_frames_in_dataset]
        if retval > lo + n_frames
            retval = 0
        end
    end

    retval
end
function get_frame_range(trajdata::TrajdataRaw, carid::Int)
    lo = trajdata.car2start[carid]
    framestart = trajdata.df[lo, :frame]

    n_frames = trajdata.df[lo, :n_frames_in_dataset]
    frameend = framestart + n_frames - 1

    framestart:frameend
end
function get_state_list(trajdata::TrajdataRaw, carid::Int, frames::AbstractVector{Int})
    veh = Vehicle()
    retval = Array(VehicleState, length(frames))
    for (i, frame) in enumerate(frames)
        get_vehicle!(veh, trajdata, carid, frame)
        retval[i] = veh.state
    end
    retval
end
function get_state_list_global(trajdata::TrajdataRaw, carid::Int, frames::AbstractVector{Int})
    veh = Vehicle()
    retval = Array(VecSE2, length(frames))
    for (i, frame) in enumerate(frames)
        get_vehicle!(veh, trajdata, carid, frame)
        retval[i] = veh.state.posG
    end
    retval
end
function get_state_list_frenet(trajdata::TrajdataRaw, carid::Int, frames::AbstractVector{Int})
    veh = Vehicle()
    retval = Array(Frenet, length(frames))
    for (i, frame) in enumerate(frames)
        get_vehicle!(veh, trajdata, carid, frame)
        retval[i] = veh.state.posF
    end
    retval
end

function pull_vehicle_headings!(trajdata::TrajdataRaw;
    v_cutoff::Float64 = 2.5, # speeds below this will use a linearly interpolated heading [fps]
    smoothing_width::Float64 = 0.5, # [s]
    )

    df = trajdata.df

    for carid in carid_set(trajdata)

        frames = collect(get_frame_range(trajdata, carid))
        states = get_state_list_global(trajdata, carid, frames)

        arr_x = map(p->p.x, states)-states[1].x
        arr_y = map(p->p.y, states)-states[1].y

        arr_v = map(i->hypot(states[i]-states[i-1]), 2:length(states)) / 0.1
        slow_indeces = find(v->v < v_cutoff, arr_v)

        slow_segments = Tuple{Int,Int}[]
        i = 0
        while i < length(slow_indeces)
            i += 1
            i_lo = i
            while i < length(slow_indeces) && slow_indeces[i+1] == slow_indeces[i]+1
                i += 1
            end
            i_hi = i
            push!(slow_segments, (slow_indeces[i_lo], slow_indeces[i_hi]))
        end
        slow_segments

        arr_heading = map(i->atan2(convert(VecE2, states[i]-states[i-1])), 2:length(states))
        unshift!(arr_heading, arr_heading[1])

        arr_x_smoothed = symmetric_exponential_moving_average(arr_x, smoothing_width)
        arr_y_smoothed = symmetric_exponential_moving_average(arr_y, smoothing_width)
        arr_dx_smoothed = arr_x_smoothed[2:end] - arr_x_smoothed[1:end-1]
        arr_dy_smoothed = arr_y_smoothed[2:end] - arr_y_smoothed[1:end-1]
        arr_heading2 = map(i->atan2(arr_y_smoothed[i]-arr_y_smoothed[i-1], arr_x_smoothed[i]-arr_x_smoothed[i-1]), 2:length(states))
        unshift!(arr_heading2, arr_heading2[1])

        arr_heading3 = deepcopy(arr_heading2)
        for (lo,hi) in slow_segments
            heading_lo = arr_heading2[max(lo-1, 1)]
            heading_hi = arr_heading2[min(hi+1, length(arr_heading2))]
            for f in lo : hi
                arr_heading3[f] = heading_lo + (heading_hi - heading_lo)*(f-lo)/(hi-lo+1)
            end
        end

        for (i,frame) in enumerate(frames)
            df_ind = car_df_index(trajdata, carid, frame)
            df[df_ind, :global_heading] = arr_heading3[i]
        end
    end

    @assert(findfirst(v->isnan(v), df[:global_heading]) == 0)

    trajdata
end
function pull_vehicle_frenet_data!(trajdata::TrajdataRaw)

    df = trajdata.df

    for i in 1 : nrow(df)

        carid = df[i, :id]
        frame = df[i, :frame]

        x = df[i, :global_x]
        y = df[i, :global_y]
        θ = df[i, :global_heading]

        posG = VecSE2(x,y,θ)
        frenet = project_posG_to_frenet!(posG, trajdata.roadway)

        df[i, :laneid] = frenet.laneid
        df[i, :frenet_extind] = frenet.extind
        df[i, :frenet_s] = frenet.s
        df[i, :frenet_t] = frenet.t
        df[i, :frenet_heading] = frenet.ϕ
    end

    trajdata
end

type FilterTrajectoryResult
    carid::Int
    x_arr::Vector{Float64}
    y_arr::Vector{Float64}
    θ_arr::Vector{Float64}
    v_arr::Vector{Float64}

    function FilterTrajectoryResult(trajdata::TrajdataRaw, carid::Int)
        dfstart = trajdata.car2start[carid]
        N = trajdata.df[dfstart, :n_frames_in_dataset]

        # these are our observations
        x_arr = fill(NaN, N)
        y_arr = fill(NaN, N)
        θ_arr = fill(NaN, N)
        v_arr = fill(NaN, N)

        for i in 1 : N
            x_arr[i] = trajdata.df[dfstart + i - 1, :global_x]
            y_arr[i] = trajdata.df[dfstart + i - 1, :global_y]
        end

        # choose an initial belief
        θ_arr[1] = atan2(y_arr[5] - y_arr[1], x_arr[5] - x_arr[1])
        v_arr[1] = trajdata.df[dfstart, :speed] #hypot(ftr.y_arr[lookahead] - y₀, ftr.x_arr[lookahead] - x₀)/ν.Δt
        if v_arr[1] < 1.0 # small speed
            # estimate with greater lookahead
            θ_arr[1] = atan2(y_arr[end] - y_arr[1], x_arr[end] - x_arr[1])
        end

        new(carid, x_arr, y_arr, θ_arr, v_arr)
    end
end
Base.length(ftr::FilterTrajectoryResult) = length(ftr.x_arr)
function filter_trajectory!(ftr::FilterTrajectoryResult, ν::VehicleSystem = VehicleSystem())

    μ = [ftr.x_arr[1], ftr.y_arr[1], ftr.θ_arr[1], ftr.v_arr[1]]

    σ = 1e-1
    Σ = diagm([σ*0.01, σ*0.01, σ*0.1, σ])

    # assume control is centered
    u = [0.0, 0.0]
    z = [NaN, NaN]

    for i in 2 : length(ftr)

        # pull observation
        z[1] = ftr.x_arr[i]
        z[2] = ftr.y_arr[i]

        # apply extended Kalman filter
        μ, Σ = EKF(ν, μ, Σ, u, z)

        # store results
        ftr.x_arr[i] = μ[1]
        ftr.y_arr[i] = μ[2]
        ftr.θ_arr[i] = μ[3]
        ftr.v_arr[i] = μ[4]
    end

    ftr
end
function Base.copy!(trajdata::TrajdataRaw, ftr::FilterTrajectoryResult)

    dfstart = trajdata.car2start[ftr.carid]
    N = trajdata.df[dfstart, :n_frames_in_dataset]

    # copy results back to trajdata
    for i in 1 : N
        trajdata.df[dfstart + i - 1, :global_x] = ftr.x_arr[i]
        trajdata.df[dfstart + i - 1, :global_y] = ftr.y_arr[i]
        trajdata.df[dfstart + i - 1, :speed]   = ftr.v_arr[i]
        trajdata.df[dfstart + i - 1, :global_heading] = ftr.θ_arr[i]
    end

    trajdata
end

function filter_trajectory!(trajdata::TrajdataRaw, carid::Int)
    #=
    Filters the given vehicle's trajectory using an Extended Kalman Filter
    =#

    ftr = FilterTrajectoryResult(trajdata, carid)

    # run pre-smoothing
    ftr.x_arr = symmetric_exponential_moving_average(ftr.x_arr, SMOOTHING_WIDTH_POS)
    ftr.y_arr = symmetric_exponential_moving_average(ftr.y_arr, SMOOTHING_WIDTH_POS)

    filter_trajectory!(ftr)

    copy!(trajdata, ftr)
    trajdata
end
function symmetric_exponential_moving_average!(trajdata::TrajdataRaw)

    for carid in carid_set(trajdata)

        dfstart = trajdata.car2start[carid]
        N = trajdata.df[dfstart, :n_frames_in_dataset]

        x_arr = Array(Float64, N)
        y_arr = Array(Float64, N)
        v_arr = Array(Float64, N)

        for i = 1 : N
            x_arr[i] = trajdata.df[dfstart + i - 1, :global_x]
            y_arr[i] = trajdata.df[dfstart + i - 1, :global_y]
            v_arr[i] = trajdata.df[dfstart + i - 1, :speed]
        end

        x_arr2 = symmetric_exponential_moving_average(x_arr, SMOOTHING_WIDTH_POS)
        y_arr2 = symmetric_exponential_moving_average(y_arr, SMOOTHING_WIDTH_POS)
        v_arr2 = symmetric_exponential_moving_average(v_arr, SMOOTHING_WIDTH_SPEED)

        for i = 1 : N
            trajdata.df[dfstart + i - 1, :global_x] = x_arr2[i]
            trajdata.df[dfstart + i - 1, :global_y] = y_arr2[i]
            trajdata.df[dfstart + i - 1, :speed]   = v_arr2[i]
        end
    end

    trajdata
end
function load_trajdata_raw(filepath::AbstractString)

    roadway = get_roadway_for_trajdata(filepath)

    print("loading from file: "); tic()
    trajdata = TrajdataRaw(filepath, roadway)
    toc()

    if splitext(filepath)[2] == ".txt" # txt is original
        # print("smoothing:         "); tic()
        # symmetric_exponential_moving_average!(trajdata)
        # toc()

        # print("global headings:   "); tic()
        # pull_vehicle_headings!(trajdata)
        # toc()

        print("filtering:         "); tic()
        for carid in carid_set(trajdata)
            filter_trajectory!(trajdata, carid)
        end
        toc()

        print("frenet data:       "); tic()
        pull_vehicle_frenet_data!(trajdata)
        toc()
    end

    trajdata
end
function input_path_to_extracted_trajdata_csv(input_path::AbstractString)
    file = splitdir(input_path)[2]
    joinpath("/media/tim/DATAPART1/PublicationData/2016_deepdrive/extracted_trajdata/", splitext(file)[1] * ".csv")
end

###############

type Trajdata
    id         :: Int                    # id assigned to this trajdata
    roadway    :: Roadway

    vehicles   :: Vector{Vehicle}        # id → Vehicle (w/ potentially uninitialized state)
    car2start  :: Vector{Int}            # id → starting index in the df
    frame2cars :: Vector{Vector{Int}}    # frame → list of carids in the scene

    frames              :: Vector{Int}          # [nrow df]
    n_frames_in_dataset :: Vector{Int}          # [nrow df]
    states              :: Vector{VehicleState} # [nrow df]

    function Trajdata(filepath::AbstractString, trajdata_id::Int=0)
        tdraw = load_trajdata_raw(filepath)
        df = tdraw.df

        id_map = Dict{Int,Int}() # old -> new

        car2start = Array(Int, length(tdraw.car2start))
        vehicles = Array(Vehicle, length(tdraw.car2start))
        for (id_old,dfind) in tdraw.car2start

            id_new = length(id_map)+1
            id_map[id_old] = id_new

            veh = Vehicle()

            veh.id = id_new
            veh.class  = df[dfind, :class ]
            veh.length = df[dfind, :length]
            veh.width  = df[dfind, :width ]

            vehicles[id_new] = veh
            car2start[id_new] = dfind
        end

        states = Array(VehicleState, nrow(df))
        for dfind in 1 : nrow(df)

            posG = VecSE2(df[dfind, :global_x], df[dfind, :global_y], df[dfind, :global_heading])
            posF = Frenet(df[dfind, :laneid], df[dfind, :frenet_extind], df[dfind, :frenet_s], df[dfind, :frenet_t], df[dfind, :frenet_heading])
            speed = df[dfind, :speed]

            states[dfind] = VehicleState(posG, posF, speed)
        end

        retval = new()
        retval.id = trajdata_id
        retval.vehicles = vehicles
        retval.car2start = car2start

        retval.frames = convert(Vector{Int}, df[:frame])
        frame_lo = minimum(retval.frames)
        Δframe = frame_lo -1
        retval.frames .-= Δframe # set low frame to 1

        frame_hi = maximum(retval.frames)
        retval.frame2cars = Array(Vector{Int}, frame_hi)
        for frame in 1:frame_hi
            carids = deepcopy(get(tdraw.frame2cars, frame+Δframe, Int[]))
            for (i,id_old) in enumerate(carids)
                carids[i] = id_map[id_old]
            end
            retval.frame2cars[frame] = carids
        end

        retval.roadway = tdraw.roadway
        retval.n_frames_in_dataset = convert(Vector{Int}, df[:n_frames_in_dataset])
        retval.states = states
        retval
    end
    function Trajdata(trajdata_id::Int, recompute::Bool=false)
        filepath = TRAJDATA_INPUT_PATHS[trajdata_id]
        if !recompute
            filepath = input_path_to_extracted_trajdata_csv(filepath)
        end

        Trajdata(filepath, trajdata_id)
    end
end

nframes(trajdata::Trajdata) = length(trajdata.frame2cars)
frame_inbounds(trajdata::Trajdata, frame::Int) = 1 ≤ frame ≤ length(trajdata.frame2cars)
carsinframe(trajdata::Trajdata, frame::Int) = trajdata.frame2cars[frame]
nth_carid(trajdata::Trajdata, frame::Int, n::Int) = trajdata.frame2cars[frame][n]
first_carid(trajdata::Trajdata, frame::Int) = nth_carid(trajdata, frame, 1)
iscarinframe(trajdata::Trajdata, carid::Int, frame::Int) = in(carid, trajdata.frame2cars[frame])

function car_df_index(trajdata::Trajdata, carid::Int, frame::Int)
    #=
    given frame and carid, find index of car in trajdata
    Returns 0 if it does not exist
    =#

    lo = trajdata.car2start[carid]
    framestart = trajdata.frames[lo]

    retval = frame - framestart + lo
    n_frames = trajdata.n_frames_in_dataset[lo]
    if retval > lo + n_frames
        retval = 0
    end

    retval
end
function get_frame_range(trajdata::Trajdata, carid::Int)
    lo = trajdata.car2start[carid]
    framestart = trajdata.frames[lo]
    n_frames = trajdata.n_frames_in_dataset[lo]
    frameend = framestart + n_frames - 1

    framestart:frameend
end
function get_vehiclestate(
    trajdata::Trajdata,
    carid::Int,
    frame::Int,
    )

    dfind = car_df_index(trajdata, carid, frame)
    trajdata.states[dfind]
end
function get_vehicle(
    trajdata::Trajdata,
    carid::Int,
    frame::Int,
    )

    veh = trajdata.vehicles[carid]
    dfind = car_df_index(trajdata, carid, frame)
    veh.state = trajdata.states[dfind]

    veh
end

function get_acceleration(trajdata::Trajdata, id::Int, frame::Int)
    if frame == 1 || !frame_inbounds(trajdata, frame)
        return 0.0 # no past info, assume zero
    end

    v_past = get_vehiclestate(trajdata, id, frame-1).v
    v_curr = get_vehiclestate(trajdata, id, frame).v

    (v_curr - v_past) / NGSIM_TIMESTEP # [ft/s²]
end
function get_turnrate(trajdata::Trajdata, id::Int, frame::Int, frenet::Bool=false)
    if frame == 1 || !frame_inbounds(trajdata, frame)
        return 0.0 # no past info, assume zero
    end

    if frenet
        past = get_vehiclestate(trajdata, id, frame-1).posF.ϕ
        curr = get_vehiclestate(trajdata, id, frame).posF.ϕ
        (curr - past) / NGSIM_TIMESTEP # [ft/s²]
    else # global frame
        past = get_vehiclestate(trajdata, id, frame-1).posG.θ
        curr = get_vehiclestate(trajdata, id, frame).posG.θ
        (curr - past) / NGSIM_TIMESTEP # [ft/s²]
    end
end