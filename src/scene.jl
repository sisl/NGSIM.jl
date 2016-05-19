type Scene
    roadway_name::Symbol
    vehicles::Vector{Vehicle} # this is a pre-allocated array that is at least as large as the maximum number of vehicles in a Trajdata frame
    n_vehicles::Int

    Scene(n_vehicles::Int=500) = new(:none, Array(Vehicle, n_vehicles), 0)
    function Scene(
        roadway_name::Symbol,
        vehicles::Vector{Vehicle},
        n_vehicles::Int=length(vehicles),
        )

        new(roadway_name, vehicles, n_vehicles)
    end
end

# iteration
Base.start(scene::Scene) = 1
Base.done(scene::Scene, i::Int) = i > length(scene)
Base.next(scene::Scene, i::Int) = (scene.vehicles[i], i+1)

# copying
function Base.copy!(dest::Scene, src::Scene)
    copy!(dest.vehicles, 1, src.vehicles, 1, src.n_vehicles)
    dest.n_vehicles = src.n_vehicles
    dest.roadway_name = src.roadway_name
    dest
end

Base.length(scene::Scene) = scene.n_vehicles
Base.getindex(scene::Scene, i::Int) = scene.vehicles[i]
function Base.setindex!(scene::Scene, veh::Vehicle, i::Int)
    scene.vehicles[i] = veh
    scene
end
function Base.empty!(scene::Scene)
    scene.n_vehicles = 0
    scene
end
function Base.get!(scene::Scene, trajdata::Trajdata, frame::Int)

    scene.roadway_name = trajdata.roadway.name

    if frame_inbounds(trajdata, frame)
        carids = trajdata.frame2cars[frame]
        scene.n_vehicles = length(carids)
        for i in 1 : scene.n_vehicles
            scene.vehicles[i] = get_vehicle(trajdata, carids[i], frame)
        end
    else
        scene.n_vehicles = 0
    end

    scene
end

get_roadway(scene::Scene) = ROADWAY_DICT[scene.roadway_name]

function get_index_of_first_vehicle_with_id(scene::Scene, id::Int)
    retval = 0
    for i in 1 : scene.n_vehicles
        if scene.vehicles[i].id == id
            retval = i
            break
        end
    end
    retval
end
function get_neighbor_index_fore(scene::Scene, vehicle_index::Int;
    max_dist = 1000.0 # [ft]
    )

    #=
    This finds the closest vehicle in the same lane based on footpoint and lanetag
    Quits once it reaches a max distance
    =#

    roadway = get_roadway(scene)

    best_index = 0
    best_dist  = Inf

    ego_state = scene.vehicles[vehicle_index].state
    active_laneid = ego_state.posF.laneid
    dist = -ego_state.posF.s # [m] dist along curve from host inertial to base of footpoint

    # walk forwards along the lanetag until we find a car in it or reach max dist
    while true

        for test_vehicle_index in 1 : scene.n_vehicles

            if test_vehicle_index == vehicle_index
                continue
            end

            veh_target = scene.vehicles[test_vehicle_index]
            laneid_target = veh_target.state.posF.laneid
            if laneid_target == active_laneid

                footpoint_s_target = veh_target.state.posF.s
                cand_dist = dist + footpoint_s_target

                if 0.0 < cand_dist < best_dist
                    best_dist = cand_dist
                    best_index = test_vehicle_index
                end
            end
        end

        break
    end

    best_index
end
function get_neighbor_index_rear(scene::Scene, vehicle_index::Int;
    max_dist = 1000.0 # [ft]
    )

    roadway = get_roadway(scene)

    best_index = 0
    best_dist   = Inf

    veh = scene.vehicles[vehicle_index]
    active_laneid = veh.state.posF.laneid
    dist = veh.state.posF.s # [m] dist along curve from host inertial to base of footpoint

    # walk forwards along the lanetag until we find a car in it or reach max dist
    while true
        for (test_vehicle_index, veh_target) in enumerate(scene)

            if test_vehicle_index == vehicle_index
                continue
            end

            laneid_target = veh_target.state.posF.laneid
            if laneid_target == active_laneid

                footpoint_s_target = veh_target.state.posF.s
                cand_dist = dist - footpoint_s_target
                if 0.0 < cand_dist < best_dist
                    best_dist = cand_dist
                    best_index = test_vehicle_index
                end
            end
        end

        break
    end

    best_index
end
function get_neighbor_index_left(scene::Scene, vehicle_index::Int, gap_lo::Float64, gap_hi::Float64,
    GAP_T::Float64 = 5.0,  # [ft]
    )

    best_index = 0
    bestΔs = Inf

    veh = scene.vehicles[vehicle_index]
    state = veh.state
    posG = state.posG
    posF = state.posF

    footpoint = get_footpoint(veh)

    right_lane_footpoint = footpoint + Vec.polar(2*GAP_T, footpoint.θ + π/2)

    for test_index in 1 : length(scene)
        if test_index != vehicle_index

            state_test = scene.vehicles[test_index].state
            if state_test.posF.laneid != posF.laneid

                # project vehicle to right lane tangent
                pos_rel = inertial2body(state_test.posG, right_lane_footpoint)
                if gap_lo ≤ pos_rel.x ≤ gap_hi && abs(pos_rel.y) ≤ GAP_T && abs(pos_rel.x) ≤ bestΔs
                    bestΔs = abs(pos_rel.x)
                    best_index = test_index
                end
            end
        end
    end

    best_index
end
function get_neighbor_index_right(scene::Scene, vehicle_index::Int, gap_lo::Float64, gap_hi::Float64,
    GAP_T::Float64 =  3*1.8,  # [ft]
    )

    #=
    Returns the index of the vehicle in the right-neighboring lane if it exists, otherwise returns 0

    - vehicle must be in the neighboring lane
    - vehicle must be within gap_lo ≤ s ≤ gap_hi
    =#

    # 1 - take the local tangent vector to your lane
    # 2 - offset it by 2d_marker_right
    # 3 - project all other vehicles to the linear-assumed lane
    # 4 - retain the vehicle that is closest in s that falls within the lane boundary

    best_index = 0
    bestΔs = Inf

    veh = scene.vehicles[vehicle_index]
    state = veh.state
    posG = state.posG
    posF = state.posF

    footpoint = get_footpoint(veh)

    right_lane_footpoint = footpoint + Vec.polar(2*GAP_T, footpoint.θ - π/2)

    for test_index in 1 : length(scene)
        if test_index != vehicle_index

            state_test = scene.vehicles[test_index].state
            if state_test.posF.laneid != posF.laneid

                # project vehicle to right lane tangent
                pos_rel = inertial2body(state_test.posG, right_lane_footpoint)
                if gap_lo ≤ pos_rel.x ≤ gap_hi && abs(pos_rel.y) ≤ GAP_T && abs(pos_rel.x) ≤ bestΔs
                    bestΔs = abs(pos_rel.x)
                    best_index = test_index
                end
            end
        end
    end

    best_index
end


function get_headway_dist_between(veh_rear::Vehicle, veh_fore::Vehicle, roadway::Roadway)

    # distance from the front of the rear vehicle to the rear of the front vehicle

    active_laneid = veh_rear.state.posF.laneid
    if veh_fore.state.posF.laneid != active_laneid
        NaN
    else
        s1 = veh_rear.state.posF.s
        s2 = veh_fore.state.posF.s
        s2 - s1 - veh_fore.length
    end
end
function get_headway_time_between(veh_rear::Vehicle, veh_fore::Vehicle, roadway::Roadway)

    active_laneid = veh_rear.state.posF.laneid

    if veh_fore.state.posF.laneid != active_laneid
        NaN
    else
        s1 = veh_rear.state.posF.s
        s2 = veh_fore.state.posF.s
        Δs = s2 - s1 - veh_fore.length
        Δs/veh_rear.state.v
    end
end
