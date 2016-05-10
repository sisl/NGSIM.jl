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

Base.start(scene::Scene) = 1
Base.done(scene::Scene, i::Int) = i > length(scene)
Base.next(scene::Scene, i::Int) = (scene.vehicles[i], i+1)

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
    curve = roadway.centerlines[active_laneid]
    dist = -(curve_at(curve, ego_state.posF.extind).s) # [m] dist along curve from host inertial to base of footpoint

    # walk forwards along the lanetag until we find a car in it or reach max dist
    while true

        for (test_vehicle_index, veh_target) in enumerate(scene)

            if test_vehicle_index == vehicle_index
                continue
            end

            laneid_target = veh_target.state.posF.laneid
            if laneid_target == active_laneid

                footpoint_s_target = curve_at(curve, veh_target.state.posF.extind).s
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
    curve = roadway.centerlines[active_laneid]
    dist = curve_at(curve, veh.state.posF.extind).s # [m] dist along curve from host inertial to base of footpoint

    # walk forwards along the lanetag until we find a car in it or reach max dist
    while true
        for (test_vehicle_index, veh_target) in enumerate(scene)

            if test_vehicle_index == vehicle_index
                continue
            end

            laneid_target = veh_target.state.posF.laneid
            if laneid_target == active_laneid

                footpoint_s_target = curve_at(curve, veh_target.state.posF.extind).s
                cand_dist = dist - footpoint_s_target
                if 0.0 < cand_dist < best_dist
                    best_dist = cand_dist
                    best_index = test_vehicle_index
                end
            end
        end

        break

        # if !isinf(best_dist) || dist > max_dist || !has_prev_lane(sn, active_lane)
        #     break
        # end

        # active_lane = prev_lane(sn, active_lane)
        # active_lanetag = active_lane.id
        # dist += active_lane.curve.s[end] # move full length
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
