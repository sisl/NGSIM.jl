const FLOATING_POINT_REGEX = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"

immutable CurvePt
    pos::VecSE2 # global position and orientation
    s::Float64 # distance along curve
end
Vec.lerp(a::CurvePt, b::CurvePt, t::Float64) = CurvePt(lerp(a.pos, b.pos, t), a.s + (b.s - a.s)*t)

immutable Frenet
    laneid::Int
    extind::Float64
    s::Float64 # distance along lane
    t::Float64 # lane offset, positive is to left
    ϕ::Float64 # lane relative heading
end
function Base.isapprox(a::Frenet, b::Frenet;
    rtol::Real=cbrt(eps(Float64)),
    atol::Real=sqrt(eps(Float64))
    )
    
    a.laneid == b.laneid &&
    isapprox(a.extind, b.extind, atol=atol, rtol=rtol) &&
    isapprox(a.s, b.s, atol=atol, rtol=rtol) &&
    isapprox(a.t, b.t, atol=atol, rtol=rtol) &&
    isapprox(a.ϕ, b.ϕ, atol=atol, rtol=rtol)
end


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
const ROADWAY_DICT = Dict{Symbol, Roadway}(ROADWAY_80.name => ROADWAY_80,
                          ROADWAY_101.name => ROADWAY_101)

function get_roadway_for_trajdata(input_path::AbstractString)
    if contains(input_path, "I-80") || contains(input_path, "trajectories-04") || contains(input_path, "trajectories-05")
        ROADWAY_80
    else
        ROADWAY_101
    end
end

function _mod2pi2(x::Float64)
    val = mod2pi(x)
    if val > pi
        val -= 2pi
    end
    return val
end
function _binary_search_curve_dist2(
    centerline::Vector{CurvePt},
    target::AbstractVec,
    )

    a = 1
    b = length(centerline)
    c = div(a+b, 2)

    @assert(length(centerline) ≥ b)

    sqdist_a = abs2(centerline[a].pos - target)
    sqdist_b = abs2(centerline[b].pos - target)
    sqdist_c = abs2(centerline[c].pos - target)

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

        right = div(c+b, 2)
        sqdist_r = abs2(centerline[right].pos - target)

        if sqdist_l < sqdist_r
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

    a # dead code
end
function _proj_rel( P₀::VecE2, P₁::VecE2, Q::VecE2 )

    #=
    Project the point (Q - P₀) onto (P₁ - P₀) and return the relative distance along
    =#

    a = Q - P₀
    b = P₁ - P₀

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

function curve_at(curve::Vector{CurvePt}, extind::Float64)
    ind_lo, ind_hi = _get_ind_lo_and_hi(curve, extind)
    curvept_lo = curve[ind_lo]
    curvept_hi = curve[ind_hi]
    t = extind - ind_lo
    lerp(curvept_lo, curvept_hi, t)
end
function get_extind(curve::Vector{CurvePt}, s::Float64)

    # get the extind for the closest s-location on the curve

    if s < 0.0
        return 1.0
    elseif s > curve[end].s
        return convert(Float64, length(curve))
    end

    a = 1
    b = length(curve)

    fa = curve[a].s - s
    fb = curve[b].s - s

    n = 1
    while true
        if b == a+1
            return a + -fa/(fb-fa)
        end

        c = div(a+b, 2)
        fc = curve[c].s - s
        n += 1

        if sign(fc) == sign(fa)
            a, fa = c, fc
        else
            b, fb = c, fc
        end
    end

    error("get_extind failed")
    NaN
end
function move_extind_along(extind::Float64, curve::Vector{CurvePt}, Δs::Float64)

    L = length(curve)
    ind_lo, ind_hi = _get_ind_lo_and_hi(curve, extind)

    s_lo = curve[ind_lo].s
    s_hi = curve[ind_hi].s
    s = lerp(s_lo, s_hi, extind-ind_lo)

    if Δs > 0.0

        if s + Δs > s_hi && ind_hi < L
            while s + Δs > s_hi && ind_hi < L
                Δs -= (s_hi - s)
                s = s_hi
                ind_lo += 1
                ind_hi += 1
                s_lo = curve[ind_lo].s
                s_hi = curve[ind_hi].s
            end
        else
            Δs = s + Δs - s_lo
        end

        t = Δs/(s_hi - s_lo)
        extind = lerp(ind_lo, ind_hi, t)
    elseif Δs < 0.0
        if s + Δs < s_lo  && ind_lo > 1
            while s + Δs < s_lo  && ind_lo > 1
                Δs += (s - s_lo)
                s = s_lo
                ind_lo -= 1
                ind_hi -= 1
                s_lo = curve[ind_lo].s
                s_hi = curve[ind_hi].s
            end
        else
            Δs = s + Δs - s_lo
        end

        t = 1.0 - Δs/(s_hi - s_lo)
        extind = lerp(ind_lo, ind_hi, t)
    end

    clamp(extind, 1.0, L)
end


immutable CurveProjection
    extind::Float64
    t::Float64 # lane offset
    ϕ::Float64 # lane-relative heading
end
function project_to_lane(posG::VecSE2, curve::Vector{CurvePt})

    # 1 - find the index of the point closest to the curve
    ind = _binary_search_curve_dist2(curve, posG)

    # 2 - interpolate between points
    extind_curve = 0.0
    p_curve = VecSE2(NaN, NaN, NaN)
    d = NaN

    if ind > 1 && ind < length(curve)
        t_lo = _proj_rel( curve[ind-1], curve[ind],   posG )
        t_hi = _proj_rel( curve[ind],   curve[ind+1], posG )

        p_lo = lerp( curve[ind-1].pos, curve[ind].pos,   t_lo )
        p_hi = lerp( curve[ind].pos,   curve[ind+1].pos, t_hi )

        d_lo = hypot( p_lo - posG )
        d_hi = hypot( p_hi - posG )

        if d_lo < d_hi
            p_curve = p_lo
            d = d_lo
            extind_curve = ind-1+t_lo
        else
            p_curve = p_hi
            d = d_hi
            extind_curve = ind+t_hi
        end
    elseif ind == 1
        t = _proj_rel( curve[1], curve[2], posG )
        p_curve = lerp( curve[1].pos, curve[2].pos, t)
        d = hypot(p_curve - posG)
        extind_curve = ind + t
    else
        t = _proj_rel( curve[end-1], curve[end], posG )
        p_curve = lerp( curve[end-1].pos, curve[end].pos, t)
        d = hypot(p_curve - posG)
        extind_curve = ind-1 + t
    end

    # 3 - compute frenet value
    dyaw = _mod2pi2( atan2( posG - p_curve ) - posG.θ )

    on_left_side = abs(_mod2pi2(dyaw - pi/2)) < abs(_mod2pi2(dyaw - 3pi/2))
    d *= on_left_side ? 1.0 : -1.0 # left side is positive, right side is negative

    yaw = _mod2pi2(posG.θ-p_curve.θ)

    CurveProjection(extind_curve, d, yaw)
end
function project_to_closest_lane(posG::VecSE2, roadway::Roadway)

    best_posF = CurveProjection(Inf, Inf, NaN)
    best_laneid = -1

    for laneid in 1 : length(roadway.centerlines)
        posF = project_to_lane(posG, roadway.centerlines[laneid])

        if abs(posF.t) < abs(best_posF.t)
            best_posF = posF
            best_laneid = laneid
        end
    end

    (best_posF, best_laneid)
end
function project_posG_to_frenet(posG::VecSE2, roadway::Roadway)
    posF, laneid = project_to_closest_lane(posG, roadway)

    extind = posF.extind
    curve = roadway.centerlines[laneid]
    s = curve_at(curve, extind).s

    Frenet(laneid, extind, s, posF.t, posF.ϕ)
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