const FLOATING_POINT_REGEX = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"

immutable CurvePt
    pos::VecSE2 # global position
    s::Float64 # distance along curve
end
Vec.lerp(a::CurvePt, b::CurvePt, t::Float64) = CurvePt(lerp(a.pos, b.pos, t), a.s + (b.s - a.s)*t)

type Roadway
    name::Symbol
    boundaries::Vector{Vector{VecE2}}
    centerlines::Vector{Vector{CurvePt}}
end
immutable RoadwayInputParams
    filepath_boundaries::ASCIIString
    filepath_centerlines::ASCIIString
end

const ROADWAY_INPUT_80 = RoadwayInputParams(Pkg.dir("NGSIM", "data", "boundaries80.txt"),
                                            Pkg.dir("NGSIM", "data", "centerlines80.txt"))
const ROADWAY_INPUT_101 = RoadwayInputParams(Pkg.dir("NGSIM", "data", "boundaries101.txt"),
                                            Pkg.dir("NGSIM", "data", "centerlines101.txt"))

function read_boundaries(io::IO)

    lines = readlines(io)
    for (i,line) in enumerate(lines)
        lines[i] = strip(line)
    end
    @assert(lines[1] == "BOUNDARIES")

    n_boundaries = parse(Int, lines[2])
    @assert(n_boundaries ≥ 0)

    retval = Array(Vector{VecE2}, n_boundaries)
    line_index = 2
    for i in 1 : n_boundaries

        @assert(lines[line_index+=1] == @sprintf("BOUNDARY %d", i))
        npts = parse(Int, lines[line_index+=1])
        line = Array(VecE2, npts)
        for j in 1 : npts
            matches = matchall(FLOATING_POINT_REGEX, lines[line_index+=1])
            x = parse(Float64, matches[1])
            y = parse(Float64, matches[2])
            line[j] = VecE2(x,y)
        end
        retval[i] = line
    end

    retval
end
function read_centerlines(io::IO)

    lines = readlines(io)
    for (i,line) in enumerate(lines)
        lines[i] = strip(line)
    end
    @assert(lines[1] == "CENTERLINES")

    n_centerlines = parse(Int, lines[2])
    @assert(n_centerlines ≥ 0)

    line_index = 2
    retval = Dict{AbstractString, Vector{CurvePt}}()
    for i in 1 : n_centerlines
        @assert(lines[line_index+=1] == "CENTERLINE")
        name = lines[line_index+=1]
        npts = parse(Int, lines[line_index+=1])
        line = Array(VecE2, npts)
        for j in 1 : npts
            matches = matchall(FLOATING_POINT_REGEX, lines[line_index+=1])
            x = parse(Float64, matches[1])
            y = parse(Float64, matches[2])
            line[j] = VecE2(x,y)
        end

        # post-process to extract heading and distance along lane
        centerline = Array(CurvePt, npts)
        let
            θ = atan2(line[2]-line[1])
            centerline[1] = CurvePt(VecSE2(line[1],θ), 0.0)
            for i in 2 : npts-1
                s = centerline[i-1].s + hypot(line[i] - line[i-1])
                θ = (atan2(line[i]-line[i-1]) + atan2(line[i+1]-line[i]))/2 # mean angle
                centerline[i] = CurvePt(VecSE2(line[i],θ), s)
            end
            s = centerline[npts-1].s + hypot(line[npts] - line[npts-1])
            θ = atan2(line[npts]-line[npts-1])
            centerline[npts] = CurvePt(VecSE2(line[npts],θ), s)
        end

        retval[name] = centerline

    end

    retval
end
function read_roadway(input_params::RoadwayInputParams)
    boundaries = open(read_boundaries, input_params.filepath_boundaries)
    centerlines = open(read_centerlines, input_params.filepath_centerlines)

    name = symbol(splitext(splitdir(input_params.filepath_boundaries)[2])[1])
    Roadway(name, boundaries, collect(values(centerlines)))
end

const ROADWAY_80 = read_roadway(ROADWAY_INPUT_80)
const ROADWAY_101 = read_roadway(ROADWAY_INPUT_101)

function get_roadway_for_trajdata(input_path::AbstractString)
    if contains(input_path, "I-80") || contains(input_path, "trajectories-04") || contains(input_path, "trajectories-05")
        ROADWAY_80
    else
        ROADWAY_101
    end
end

function _mod2pi2(X::Float64)
    val = mod2pi(X)
    if val > pi
        val -= 2pi
    end
    return val
end
function _binary_search_curve_dist2(
    centerline::Vector{CurvePt},
    target::AbstractVec,
    sq_dist_threshold::Float64 = 0.1 # if dist is less than this we know we have the closest point
    )

    a = 1
    b = length(centerline)

    @assert(length(centerline) ≥ b)

    sqdist_a = abs2(centerline[a].pos - target)
    sqdist_b = abs2(centerline[b].pos - target)

    if b == a
        return a
    elseif b == a + 1
        return sqdist_b < sqdist_a ? b : a
    end

    c = div(a+b, 2)
    sqdist_c = abs2(centerline[c].pos - target)

    n = 1
    while true
        if b == a
            return a
        elseif b == a + 1
            return sqdist_b < sqdist_a ? b : a
        elseif c == a + 1 && c == b - 1
            if sqdist_a < sqdist_b && sqdist_a < sqdist_c
                return a
            elseif sqdist_b < sqdist_a && sqdist_b < sqdist_c
                return b
            else
                return c
            end
        end

        left = div(a+c, 2)
        sqdist_l = abs2(centerline[left].pos - target)

        if sqdist_l < sq_dist_threshold
            return left
        end

        right = div(c+b, 2)
        sqdist_r = abs2(centerline[right].pos - target)

        if sqdist_r < sq_dist_threshold
            return right
        elseif sqdist_l < sqdist_r
            b = c
            sqdist_b = sqdist_c
            c = left
            sqdist_c = sqdist_l
        else
            a = c
            sqdist_a = sqdist_c
            c = right
            sqdist_c = sqdist_r
        end
    end

    a
end

function _proj_rel( P₀::VecE2, P₁::VecE2, Q::VecE2 )

    #=
    Project the point (Q - P₀) onto (P₁ - P₀) and return the relative distance along
    =#

    b = P₁ - P₀
    a = Q - P₀

    c = proj(a, b, VecE2)

    if b.x != 0.0
        t = c.x / b.x
    else
        t = c.y / b.y
    end

    clamp(t, 0.0, 1.0)
end
_proj_rel( P₀::CurvePt, P₁::CurvePt, Q::VecSE2 ) = _proj_rel(convert(VecE2, P₀.pos), convert(VecE2, P₁.pos), convert(VecE2, Q))

function _get_ind_lo_and_hi(curve::Vector{CurvePt}, extind::Float64)
    ind_lo = floor(Int, extind)
    ind_hi = ceil(Int, extind)
    L = length(curve)

    if ind_hi < 2
        ind_lo = 1
        ind_hi = 2
    elseif ind_lo ≥ L
        ind_hi = L
        ind_lo = ind_hi - 1
    end

    if ind_lo == ind_hi
        if ind_lo == 1
            ind_hi += 1
        else
            ind_lo -= 1
        end
    end

    (ind_lo, ind_hi)
end

function project_to_lane(posG::VecSE2, curve::Vector{CurvePt})

    # 1 - find the index of the point closest to the curve
    ind = _binary_search_curve_dist2(curve, posG)

    # 2 - interpolate between points
    extind_curve = 0.0
    p_curve = VecSE2(NaN, NaN, NaN)
    if ind > 1 && ind < length(curve)
        t_lo = _proj_rel( curve[ind-1], curve[ind],   posG )
        t_hi = _proj_rel( curve[ind],   curve[ind+1], posG )

        p_lo = lerp( curve[ind-1].pos, curve[ind].pos,   t_lo )
        p_hi = lerp( curve[ind].pos,   curve[ind+1].pos, t_hi )

        d_lo = hypot( p_lo - posG )
        d_hi = hypot( p_hi - posG )

        if d_lo < d_hi
            p_curve = p_lo
            extind_curve = ind-1.0+t_lo
        else
            p_curve = p_hi
            extind_curve = ind+t_hi
        end
    elseif ind == 1
        t = _proj_rel( curve[ind], curve[ind+1], posG )
        p_curve = lerp( curve[ind].pos, curve[ind+1].pos, t)
        extind_curve = ind + t
    else
        t = _proj_rel( curve[ind-1], curve[ind], posG )
        p_curve = lerp( curve[ind-1].pos, curve[ind].pos, t)
        extind_curve = ind - 1.0 + t
    end

    # 3 - compute frenet value
    d = hypot(p_curve - posG)

    dyaw = mod2pi( atan2( posG - p_curve ) - posG.θ )

    on_left_side = abs(_mod2pi2(dyaw - pi/2)) < abs(_mod2pi2(dyaw - 3pi/2))
    d *= on_left_side ? 1.0 : -1.0 # left side is positive, right side is negative

    yaw = _mod2pi2(posG.θ-p_curve.θ)

    # 4 - return Frenet pt {extind, d, yaw}
    VecSE2(extind_curve, d, yaw)
end
function project_to_closest_lane(posG::VecSE2, roadway::Roadway)

    best_posF = VecSE2(Inf, Inf, NaN)
    best_laneid = -1

    for (laneid, curve) in enumerate(roadway.centerlines)
        posF = project_to_lane(posG, curve)

        if abs(posF.y) < abs(best_posF.y)
            best_posF = posF
            best_laneid = laneid
        end
    end

    (best_posF, best_laneid)
end
function curve_at(curve::Vector{CurvePt}, extind::Float64)
    ind_lo, ind_hi = _get_ind_lo_and_hi(curve, extind)
    curvept_lo = curve[ind_lo]
    curvept_hi = curve[ind_hi]
    t = extind - ind_lo
    lerp(curvept_lo, curvept_hi, t)
end

function move_extind_along(extind::Float64, curve::Vector{CurvePt}, Δs::Float64)

    L = length(curve)
    ind_lo, ind_hi = _get_ind_lo_and_hi(curve, extind)

    s_lo = curve[ind_lo].s
    s_hi = curve[ind_hi].s
    s = lerp(s_lo, s_hi, extind-ind_lo)

    if Δs > 0.0

        while s + Δs > s_hi && ind_hi < L
            Δs -= (s_hi - s)
            s = s_hi
            ind_lo += 1
            ind_hi += 1
            s_lo = curve[ind_lo].s
            s_hi = curve[ind_hi].s
        end

        t = Δs/(s_hi - s_lo)
        extind = lerp(ind_lo, ind_hi, t)
    elseif Δs < 0.0
        while s + Δs < s_lo  && ind_lo > 1
            Δs += (s - s_lo)
            s = s_lo
            ind_lo -= 1
            ind_hi -= 1
            s_lo = curve[ind_lo].s
            s_hi = curve[ind_hi].s
        end

        t = 1.0 - Δs/(s_hi - s_lo)
        extind = lerp(ind_lo, ind_hi, t)
    end

    clamp(extind, 1.0, L)
end

function get_neighbor_laneid_left(roadway::Roadway, laneid::Int, extind::Float64)

    # Returns 0 on failure

    curve = roadway.centerlines[laneid]
    footpoint = curve_at(curve, extind)

    # project perpendicular to the left
    left_lane_candidate = footpoint.pos + Vec.polar(8.0, footpoint.pos.θ + π/2)

    best_posF, best_laneid = project_to_closest_lane(left_lane_candidate, roadway)

    if best_laneid == laneid
        best_laneid = 0
    end
    best_laneid
end
function get_neighbor_laneid_right(roadway::Roadway, laneid::Int, extind::Float64)

    # Returns 0 on failure

    curve = roadway.centerlines[laneid]
    footpoint = curve_at(curve, extind)

    # project perpendicular to the right
    right_lane_candidate = footpoint.pos + Vec.polar(8.0, footpoint.pos.θ - π/2)
    best_posF, best_laneid = project_to_closest_lane(right_lane_candidate, roadway)

    if best_laneid == laneid
        best_laneid = 0
    end
    best_laneid
end