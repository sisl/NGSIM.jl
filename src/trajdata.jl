const NGSIM_TIMESTEP = 0.1 # [sec]
const SMOOTHING_WIDTH_POS = 0.5 # [s]

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

###############

type NGSIMTrajdata
    df         :: DataFrame
    car2start  :: Dict{Int, Int}         # maps carindex to starting index in the df
    frame2cars :: Dict{Int, Vector{Int}} # maps frame to list of carids in the scene

    function NGSIMTrajdata(input_path::AbstractString)

        @assert(isfile(input_path))

        df = readtable(input_path, separator=' ', header = false)
        col_names = [:id, :frame, :n_frames_in_dataset, :epoch, :local_x, :local_y, :global_x, :global_y, :length, :width, :class, :speed, :acc, :lane, :carind_front, :carind_rear, :dist_headway, :time_headway]
        for (i,name) in enumerate(col_names)
            rename!(df, Symbol(@sprintf("x%d", i)), name)
        end

        df[:global_heading] = fill(NaN, nrow(df))

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

        new(df, car2start, frame2cars)
    end
end

Records.nframes(trajdata::NGSIMTrajdata) = maximum(keys(trajdata.frame2cars))
carsinframe(trajdata::NGSIMTrajdata, frame::Int) = get(trajdata.frame2cars, frame, Int[]) # NOTE: memory allocation!
carid_set(trajdata::NGSIMTrajdata) = Set(keys(trajdata.car2start)) # NOTE: memory allocation!
nth_carid(trajdata::NGSIMTrajdata, frame::Int, n::Int) = trajdata.frame2cars[frame][n]
first_carid(trajdata::NGSIMTrajdata, frame::Int) = nth_carid(trajdata, frame, 1)
iscarinframe(trajdata::NGSIMTrajdata, carid::Int, frame::Int) = in(carid, carsinframe(trajdata, frame))

function car_df_index(trajdata::NGSIMTrajdata, carid::Int, frame::Int)
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
function get_frame_range(trajdata::NGSIMTrajdata, carid::Int)
    lo = trajdata.car2start[carid]
    framestart = trajdata.df[lo, :frame]

    n_frames = trajdata.df[lo, :n_frames_in_dataset]
    frameend = framestart + n_frames - 1

    framestart:frameend
end

function pull_vehicle_headings!(trajdata::NGSIMTrajdata;
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

type FilterTrajectoryResult
    carid::Int
    x_arr::Vector{Float64}
    y_arr::Vector{Float64}
    θ_arr::Vector{Float64}
    v_arr::Vector{Float64}

    function FilterTrajectoryResult(trajdata::NGSIMTrajdata, carid::Int)
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
function Base.copy!(trajdata::NGSIMTrajdata, ftr::FilterTrajectoryResult)

    dfstart = trajdata.car2start[ftr.carid]
    N = trajdata.df[dfstart, :n_frames_in_dataset]

    # copy results back to trajdata
    for i in 1 : N
        trajdata.df[dfstart + i - 1, :global_x] = ftr.x_arr[i]
        trajdata.df[dfstart + i - 1, :global_y] = ftr.y_arr[i]
        # trajdata.df[dfstart + i - 1, :speed]   = ftr.v_arr[i]
        if i > 1
            trajdata.df[dfstart + i - 1, :speed]   = hypot(ftr.x_arr[i] - ftr.x_arr[i-1], ftr.y_arr[i] - ftr.y_arr[i-1]) / NGSIM_TIMESTEP
        else
            trajdata.df[dfstart + i - 1, :speed]   = hypot(ftr.x_arr[i+1] - ftr.x_arr[i], ftr.y_arr[i+1] - ftr.y_arr[i]) / NGSIM_TIMESTEP
        end
        trajdata.df[dfstart + i - 1, :global_heading] = ftr.θ_arr[i]
    end

    trajdata
end

function filter_trajectory!(trajdata::NGSIMTrajdata, carid::Int)
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
function load_ngsim_trajdata(filepath::AbstractString)

    print("loading from file: "); tic()
    tdraw = NGSIMTrajdata(filepath)
    toc()

    if splitext(filepath)[2] == ".txt" # txt is original
        print("filtering:         "); tic()
        for carid in carid_set(tdraw)
            filter_trajectory!(tdraw, carid)
        end
        toc()
    end

    tdraw
end

function Base.convert(::Type{Trajdata}, tdraw::NGSIMTrajdata, roadway::Roadway)

    df = tdraw.df

    vehdefs = Dict{Int, VehicleDef}()
    states = Array(RecordState{VehicleState, Int}, nrow(df))
    frames = Array(RecordFrame, nframes(tdraw))

    for (id, dfind) in tdraw.car2start
        vehdefs[id] = VehicleDef(df[dfind, :class], df[dfind, :length]*METERS_PER_FOOT, df[dfind, :width]*METERS_PER_FOOT)
    end

    state_ind = 0
    for frame in 1 : nframes(tdraw)

        frame_lo = state_ind+1

        for id in carsinframe(tdraw, frame)
            dfind = car_df_index(tdraw, id, frame)

            posG = VecSE2(df[dfind, :global_x]*METERS_PER_FOOT, df[dfind, :global_y]*METERS_PER_FOOT, df[dfind, :global_heading])
            speed = df[dfind, :speed]*METERS_PER_FOOT

            states[state_ind += 1] = RecordState(VehicleState(posG, roadway, speed), id)
        end

        frame_hi = state_ind
        frames[frame] = RecordFrame(frame_lo, frame_hi)
    end

    Trajdata(NGSIM_TIMESTEP, frames, states, vehdefs)
end

get_corresponding_roadway(filename::String) = contains(filename, "i101") ? ROADWAY_101 : ROADWAY_80


function convert_raw_ngsim_to_trajdatas()
    for filename in ("i101_trajectories-0750am-0805am.txt",
                     "i101_trajectories-0805am-0820am.txt",
                     "i101_trajectories-0820am-0835am.txt",
                     "i80_trajectories-0400-0415.txt",
                     "i80_trajectories-0500-0515.txt",
                     "i80_trajectories-0515-0530.txt")

        println("converting ", filename); tic()

        filepath = Pkg.dir("NGSIM", "data", filename)
        roadway = get_corresponding_roadway(filename)
        tdraw = NGSIM.load_ngsim_trajdata(filepath)
        trajdata = convert(Trajdata, tdraw, roadway)

        outpath = Pkg.dir("NGSIM", "data", "trajdata_"*filename)
        open(io->write(io, MIME"text/plain"(), trajdata), outpath, "w")

        toc()
    end
end

const TRAJDATA_PATHS = [
                        Pkg.dir("NGSIM", "data", "trajdata_i101_trajectories-0750am-0805am.txt"),
                        Pkg.dir("NGSIM", "data", "trajdata_i101_trajectories-0805am-0820am.txt"),
                        Pkg.dir("NGSIM", "data", "trajdata_i101_trajectories-0820am-0835am.txt"),
                        Pkg.dir("NGSIM", "data", "trajdata_i80_trajectories-0400-0415.txt"),
                        Pkg.dir("NGSIM", "data", "trajdata_i80_trajectories-0500-0515.txt"),
                        Pkg.dir("NGSIM", "data", "trajdata_i80_trajectories-0515-0530.txt"),
                       ]

function load_trajdata(filepath::AbstractString)
    td = open(io->read(io, MIME"text/plain"(), Trajdata), filepath, "r")
    td
end
load_trajdata(i::Int) = load_trajdata(TRAJDATA_PATHS[i])
get_corresponding_roadway(i::Int) = get_corresponding_roadway(TRAJDATA_PATHS[i])
