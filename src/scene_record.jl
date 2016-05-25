type SceneRecord
    roadway_name::Symbol
    vehicles::Dict{Int, Vehicle} # vehicle definitions (id -> Vehicle)
    n_vehicles::Vector{Int} # number of vehicles in each column of states
    ids::Matrix{Int} # vehicle id in each slot
    states::Matrix{VehicleState} # vehicles × history
                                 # states[i,j] is the ith vehicle in the jth most recent frame

    function SceneRecord(history::Int, max_n_vehicles::Int=500, roadway_name::Symbol=:none)
        retval = new()
        retval.roadway_name = roadway_name
        retval.vehicles = Dict{Int, Vehicle}()
        retval.n_vehicles = zeros(Int, history)
        retval.ids = Array(Int, max_n_vehicles, history)
        retval.states = Array(VehicleState, max_n_vehicles, history)
        retval
    end
end

record_length(rec::SceneRecord) = size(rec.states, 2)
_state_ind(pastframe::Int) = 1 - pastframe

function Base.empty!(rec::SceneRecord)
    fill!(rec.n_vehicles, 0)
    rec
end

function Base.getindex(rec::SceneRecord, i::Int, pastframe::Int)
    j = _state_ind(pastframe)
    id = rec.ids[i,j]
    veh = rec.vehicles[id]
    veh.state = rec.states[i,j]
    veh
end
Base.getindex(rec::SceneRecord, i::Int) = Base.getindex(rec, i, 0)

function Base.setindex!(rec::SceneRecord, veh::Vehicle, i::Int, pastframe::Int)

    if !haskey(rec.vehicles, veh.id)
        rec.vehicles[veh.id] = deepcopy(veh)
    end

    j = _state_ind(pastframe)
    rec.ids[i,j] = veh.id
    rec.states[i,j] = veh.state

    rec
end
Base.setindex!(rec::SceneRecord, veh::Vehicle, i::Int) = Base.setindex!(rec, veh, i, 0)

get_roadway(rec::SceneRecord) = ROADWAY_DICT[rec.roadway_name]
function iscarinframe(rec::SceneRecord, id::Int, pastframe::Int=0)
    j = _state_ind(pastframe)
    for i in 1 : rec.n_vehicles[j]
        if rec.ids[i,j] == id
            return true
        end
    end
    false
end
function get_index_of_first_vehicle_with_id(rec::SceneRecord, id::Int, pastframe::Int=0)
    j = _state_ind(pastframe)
    for i in 1 : rec.n_vehicles[j]
        if rec.ids[i,j] == id
            return i
        end
    end
    0
end
function get_vehicle(rec::SceneRecord, id::Int, pastframe::Int=0)
    j = _state_ind(pastframe)
    for i in 1 : rec.n_vehicles[j]
        if rec.ids[i,j] == id
            veh = rec.vehicles[id]
            veh.state = rec.states[i,j]
            return veh
        end
    end
    error("NGSIM.get_vehicle(rec::SceneRecord) vehicle not found")
end
function get_vehiclestate(rec::SceneRecord, id::Int, pastframe::Int=0)
    j = _state_ind(pastframe)
    for i in 1 : rec.n_vehicles[j]
        if rec.ids[i,j] == id
            return rec.states[i,j]
        end
    end
    error("vehiclestate not found")
end

function push_back_records!(rec::SceneRecord)
    for i in size(rec.states, 2) : -1 : 2
        for j in 1 : rec.n_vehicles[i-1]
            rec.states[j,i] = rec.states[j,i-1]
            rec.ids[j,i] = rec.ids[j,i-1]
        end
        rec.n_vehicles[i] = rec.n_vehicles[i-1]
    end
    rec.n_vehicles[1] = 0
    rec
end
function Base.insert!(rec::SceneRecord, scene::Scene, pastframe::Int=0)
    rec.roadway_name = scene.roadway_name

    j = _state_ind(pastframe)
    for (i,veh) in enumerate(scene)
        if !haskey(rec.vehicles, veh.id)
            rec.vehicles[veh.id] = Vehicle(veh.id, veh.class, veh.length, veh.width)
        end

        rec.states[i,j] = veh.state
        rec.ids[i,j] = veh.id
    end
    rec.n_vehicles[j] = length(scene)

    rec
end
function update!(rec::SceneRecord, scene::Scene)
    push_back_records!(rec)
    insert!(rec, scene, 0)
    rec
end

function Base.get!(scene::Scene, rec::SceneRecord, pastframe::Int=0)
    scene.roadway_name = rec.roadway_name

    j = _state_ind(pastframe)
    for i in 1 : rec.n_vehicles[j]
        id = rec.ids[i,j]
        veh = deepcopy(rec.vehicles[id])
        veh.state = rec.states[i,j]
        scene.vehicles[i] = veh
    end

    scene.n_vehicles = rec.n_vehicles[j]

    scene
end

function get_neighbor_index_fore(rec::SceneRecord, vehicle_index::Int, pastframe::Int=0)

    #=
    This finds the closest vehicle in the same lane based on footpoint and lanetag
    =#

    roadway = get_roadway(rec)

    j = _state_ind(pastframe)
    ego_state = rec.states[vehicle_index, j]
    active_laneid = ego_state.posF.laneid
    dist = -ego_state.posF.s # [m] dist along curve from host inertial to base of footpoint

    # walk forwards along the lanetag until we find a car in it or reach max dist
    # TODO: link to next centerline

    best_index = 0
    best_dist  = Inf

    for test_index in 1 : rec.n_vehicles[j]

        if test_index == vehicle_index
            continue
        end

        posF_target = rec.states[test_index,j].posF
        if posF_target.laneid == active_laneid

            cand_dist = dist + posF_target.s

            if 0.0 < cand_dist < best_dist
                best_dist = cand_dist
                best_index = test_index
            end
        end
    end

    best_index
end
function get_neighbor_index_rear(rec::SceneRecord, vehicle_index::Int, pastframe::Int=0)

    roadway = get_roadway(rec)

    j = _state_ind(pastframe)
    ego_state = rec.states[vehicle_index, j]
    active_laneid = ego_state.posF.laneid
    dist = ego_state.posF.s # [m] dist along curve from host inertial to base of footpoint

    best_index = 0
    best_dist  = Inf

    for test_index in 1 : rec.n_vehicles[j]

        if test_index == vehicle_index
            continue
        end

        posF_target = rec.states[test_index,j].posF
        if posF_target.laneid == active_laneid

            cand_dist = dist - posF_target.s

            if 0.0 < cand_dist < best_dist
                best_dist = cand_dist
                best_index = test_index
            end
        end
    end

    best_index
end
function get_neighbor_index_left(rec::SceneRecord, vehicle_index::Int, gap_lo::Float64, gap_hi::Float64, pastframe::Int=0;
    GAP_T::Float64 = 5.0,  # [ft]
    )

    best_index = 0
    bestΔs = Inf

    j = _state_ind(pastframe)
    id = rec.ids[vehicle_index, j]
    veh = rec.vehicles[id]
    veh.state = rec.states[vehicle_index, j]
    posF = veh.state.posF

    footpoint = get_footpoint(veh)
    right_lane_footpoint = footpoint + Vec.polar(2*GAP_T, footpoint.θ + π/2)

    for test_index in 1 : rec.n_vehicles[j]
        if test_index == vehicle_index
            continue
        end

        state_test = rec.states[test_index,j]
        if state_test.posF.laneid != posF.laneid

            # project vehicle to right lane tangent
            pos_rel = inertial2body(state_test.posG, right_lane_footpoint)
            if gap_lo ≤ pos_rel.x ≤ gap_hi && abs(pos_rel.y) ≤ GAP_T && abs(pos_rel.x) ≤ bestΔs
                bestΔs = abs(pos_rel.x)
                best_index = test_index
            end
        end
    end

    best_index
end
function get_neighbor_index_right(rec::SceneRecord, vehicle_index::Int, gap_lo::Float64, gap_hi::Float64, pastframe::Int=0;
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

    j = _state_ind(pastframe)
    id = rec.ids[vehicle_index, j]
    veh = rec.vehicles[id]
    veh.state = rec.states[vehicle_index, j]
    posF = veh.state.posF

    footpoint = get_footpoint(veh)
    right_lane_footpoint = footpoint + Vec.polar(2*GAP_T, footpoint.θ - π/2)

    for test_index in 1 : rec.n_vehicles[j]
        if test_index == vehicle_index
            continue
        end

        state_test = rec.states[test_index,j]
        if state_test.posF.laneid != posF.laneid

            # project vehicle to right lane tangent
            pos_rel = inertial2body(state_test.posG, right_lane_footpoint)
            if gap_lo ≤ pos_rel.x ≤ gap_hi && abs(pos_rel.y) ≤ GAP_T && abs(pos_rel.x) ≤ bestΔs
                bestΔs = abs(pos_rel.x)
                best_index = test_index
            end
        end
    end

    best_index
end

function get_acceleration(rec::SceneRecord, id::Int, pastframe::Int=0)
    if pastframe > record_length(rec) || pastframe > 0 || !iscarinframe(rec, id, pastframe-1)
        return 0.0 # no past info, assume zero
    end

    v_past = get_vehiclestate(rec, id, pastframe-1).v
    v_curr = get_vehiclestate(rec, id, pastframe).v

    (v_curr - v_past) / NGSIM_TIMESTEP # [ft/s²]
end
function get_turnrate(rec::SceneRecord, id::Int, pastframe::Int=0, frenet::Bool=false)
    if pastframe > record_length(rec) || pastframe > 0 || !iscarinframe(rec, id, pastframe-1)
        return 0.0 # no past info, assume zero
    end

    if frenet
        past = get_vehiclestate(rec, id, pastframe-1).posF.ϕ
        curr = get_vehiclestate(rec, id, pastframe).posF.ϕ
        (curr - past) / NGSIM_TIMESTEP # [ft/s²]
    else # global frame
        past = get_vehiclestate(rec, id, pastframe-1).posG.θ
        curr = get_vehiclestate(rec, id, pastframe).posG.θ
        (curr - past) / NGSIM_TIMESTEP # [ft/s²]
    end
end