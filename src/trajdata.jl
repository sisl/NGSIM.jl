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
    retval = Array{Float64}(N)

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

mutable struct FilterTrajectoryResult
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
function load_ngsim_trajdata(filepath::String; autofilter::Bool=true)

    print("loading from file: "); tic()
    tdraw = NGSIMTrajdata(filepath)
    toc()

    if autofilter && splitext(filepath)[2] == ".txt" # txt is original
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
    states = Array{RecordState{VehicleState, Int}}(nrow(df))
    frames = Array{RecordFrame}(nframes(tdraw))

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
    for filename in NGSIM_TRAJDATA_PATHS
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

function load_trajdata(filepath::String)
    td = open(io->read(io, MIME"text/plain"(), Trajdata), filepath, "r")
    td
end
load_trajdata(i::Int) = load_trajdata(TRAJDATA_PATHS[i])
get_corresponding_roadway(i::Int) = get_corresponding_roadway(TRAJDATA_PATHS[i])
