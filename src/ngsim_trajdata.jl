"""
NGSIMTrajdata

The trajectory data stored in the original NGSIM dataset format.
The dataset is a white-space separated text file with columns:
    id                  - I - Vehicle identification number (ascending by time of entry into section)
    frame               - I - Frame Identification number (ascending by start time), units 1/10 of a second
    n_frames_in_dataset - I - Total number of frames in which the vehicle appears in this data set, units 1/10 of a second
    epoch               - I - Elapsed time since Jan 1, 1970, in milliseconds
    local_x             - F - Lateral (X) coordinate of the front center of the vehicle with respect to the left-most edge of the section in the direction of travel, in feet
    local_y             - F - Longitudinal (Y) coordinate of the front center of the vehicle with respect to the entry edge of the section in the direction of travel, in feet
    global_x            - F - X Coordinate of the front center of the vehicle based on CA State Plane III in NAD83
    global_y            - F - Y Coordinate of the front center of the vehicle based on CA State Plane III in NAD83
    length              - F - Length of the vehicle, in feet
    width               - F - Width of the vehicle, in feet
    class               - I - vehicle class, 1 - motorcycle, 2 - auto, 3 - truck
    speed               - F - Instantaneous velocity of vehicle, in ft/second
    acc                 - F - Instantaneous acceleration of vehicle, in ft/second^2
    lane                - I - Current lane position of vehicle
    carind_front        - I - Vehicle Id of the lead vehicle in the same lane. A value of '0' represents no preceding vehicle
    carind_rear         - I - Vehicle Id of the vehicle following the subject vehicle in the same lane. A value of '0' represents no following vehicle
    dist_headway        - F - Spacing provides the distance between the front-center of a vehicle to the front-center of the preceding vehicle, in feet
    time_headway        - F - Headway provides the time to travel from the front-center of a vehicle (at the speed of the vehicle) to the front-center of the preceding vehicle. A headway value of 9999.99 means that the vehicle is traveling at zero speed (congested conditions), in second
"""
mutable struct NGSIMTrajdata
    df         :: DataFrame
    car2start  :: Dict{Int, Int}         # maps carindex to starting index in the df
    frame2cars :: Dict{Int, Vector{Int}} # maps frame to list of carids in the scene

    function NGSIMTrajdata(input_path::String)

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

        arr_heading = map(i->atan(convert(VecE2, states[i]-states[i-1])), 2:length(states))
        unshift!(arr_heading, arr_heading[1])

        arr_x_smoothed = symmetric_exponential_moving_average(arr_x, smoothing_width)
        arr_y_smoothed = symmetric_exponential_moving_average(arr_y, smoothing_width)
        arr_dx_smoothed = arr_x_smoothed[2:end] - arr_x_smoothed[1:end-1]
        arr_dy_smoothed = arr_y_smoothed[2:end] - arr_y_smoothed[1:end-1]
        arr_heading2 = map(i->atan(arr_y_smoothed[i]-arr_y_smoothed[i-1], arr_x_smoothed[i]-arr_x_smoothed[i-1]), 2:length(states))
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

const NGSIM_TRAJDATA_PATHS = [
                        joinpath(@__DIR__, "../data/i101_trajectories-0750am-0805am.txt"),
                        joinpath(@__DIR__, "../data/i101_trajectories-0805am-0820am.txt"),
                        joinpath(@__DIR__, "../data/i101_trajectories-0820am-0835am.txt"),
                        joinpath(@__DIR__, "../data/i80_trajectories-0400-0415.txt"),
                        joinpath(@__DIR__, "../data/i80_trajectories-0500-0515.txt"),
                        joinpath(@__DIR__, "../data/i80_trajectories-0515-0530.txt"),
                       ]

load_ngsim_trajdata(i::Int) = NGSIMTrajdata(NGSIM_TRAJDATA_PATHS[i])