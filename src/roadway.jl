const FLOATING_POINT_REGEX = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
const METERS_PER_FOOT = 0.3048

mutable struct NGSIMRoadway
    name::Symbol
    boundaries::Vector{Vector{VecE2}}
    centerlines::Vector{Vector{CurvePt}}
end
struct RoadwayInputParams
    filepath_boundaries::String
    filepath_centerlines::String
end

function read_boundaries(io::IO)

    lines = readlines(io)
    for (i,line) in enumerate(lines)
        lines[i] = strip(line)
    end
    @assert(lines[1] == "BOUNDARIES")

    n_boundaries = parse(Int, lines[2])
    @assert(n_boundaries ≥ 0)

    retval = Array{Vector{VecE2}}(n_boundaries)
    line_index = 2
    for i in 1 : n_boundaries

        @assert(lines[line_index+=1] == @sprintf("BOUNDARY %d", i))
        npts = parse(Int, lines[line_index+=1])
        line = Array{VecE2}(npts)
        for j in 1 : npts
            matches = matchall(FLOATING_POINT_REGEX, lines[line_index+=1])
            x = parse(Float64, matches[1]) * METERS_PER_FOOT # convert to meters
            y = parse(Float64, matches[2]) * METERS_PER_FOOT
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
        line = Array{VecE2}(npts)
        for j in 1 : npts
            matches = matchall(FLOATING_POINT_REGEX, lines[line_index+=1])
            x = parse(Float64, matches[1]) * METERS_PER_FOOT # convert to meters
            y = parse(Float64, matches[2]) * METERS_PER_FOOT # convert to meters
            line[j] = VecE2(x,y)
        end

        # post-process to extract heading and distance along lane
        centerline = Array{CurvePt}(npts)
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
    NGSIMRoadway(name, boundaries, collect(values(centerlines)))
end

function write_lwpolyline(io::IO, pts::Vector{CurvePt}, handle_int::Int, is_closed::Bool=false)
    N = length(pts)
    println(io, "  0")
    println(io, "LWPOLYLINE")
    println(io, "  5") # handle (increases)
    @printf(io, "B%d\n", handle_int)
    println(io, "100") # subclass marker
    println(io, "AcDbEntity")
    println(io, "  8") # layer name
    println(io, "150")
    println(io, "  6") # linetype name
    println(io, "ByLayer")
    println(io, " 62") # color number
    println(io, "  256")
    println(io, "370") # lineweight enum
    println(io, "   -1")
    println(io, "100") # subclass marker
    println(io, "AcDbPolyline")
    println(io, " 90") # number of vertices
    @printf(io, "   %d\n", N)
    println(io, " 70") # 0 is default, 1 is closed
    @printf(io, "    %d\n", is_closed ? 1 : 0)
    println(io, " 43") # 0 is constant width
    println(io, "0")

    for i in 1 : N
        println(io, " 10")
        @printf(io, "%.3f\n", pts[i].pos.x)
        println(io, " 20")
        @printf(io, "%.3f\n", pts[i].pos.y)
    end
end

function convert_curves_feet_to_meters!(roadway::Roadway)
    for seg in roadway.segments
        for lane in seg.lanes
            for (i,P) in enumerate(lane.curve)
                lane.curve[i] = CurvePt(
                        VecSE2(P.pos.x*METERS_PER_FOOT, P.pos.y*METERS_PER_FOOT, P.pos.θ),
                        P.s*METERS_PER_FOOT, P.k/METERS_PER_FOOT, P.kd/METERS_PER_FOOT)
            end
        end
    end
    roadway
end

function _fit_curve(
    pts::AbstractMatrix{Float64}, # 2×n
    desired_distance_between_samples::Real;
    max_iterations::Int=50,
    epsilon::Float64=1e-4,
    n_intervals_in_arclen::Int=100,
    )

    @assert(size(pts, 1) == 2)

    spline_coeffs = fit_cubic_spline(pts)

    L = calc_curve_length(spline_coeffs[1], spline_coeffs[2], n_intervals_per_segment=n_intervals_in_arclen)
    n = round(Int, L/desired_distance_between_samples, RoundNearestTiesUp)+1

    s_arr = collect(linspace(0.0,L,n))
    t_arr = calc_curve_param_given_arclen(spline_coeffs[1], spline_coeffs[2], s_arr,
        curve_length=L, max_iterations=max_iterations, epsilon=epsilon, n_intervals_in_arclen=n_intervals_in_arclen)

    x_arr = sample_spline(spline_coeffs[1], t_arr)
    y_arr = sample_spline(spline_coeffs[2], t_arr)
    θ_arr = sample_spline_theta(spline_coeffs[1], spline_coeffs[2], t_arr)

    κ_arr = sample_spline_curvature(spline_coeffs[1], spline_coeffs[2], t_arr)
    κd_arr = sample_spline_derivative_of_curvature(spline_coeffs[1], spline_coeffs[2], t_arr)

    @assert(!any(s->isnan(s), s_arr))
    @assert(!any(s->isnan(s), x_arr))
    @assert(!any(s->isnan(s), y_arr))
    @assert(!any(s->isnan(s), θ_arr))

    curve = Array{CurvePt}(n)
    for i in 1 : n
        pos = VecSE2(x_arr[i], y_arr[i], θ_arr[i])
        curve[i] = CurvePt(pos, s_arr[i], κ_arr[i], κd_arr[i])
    end
    curve
end

"""
    read_dxf(io::IO, ::Type{Roadway})
Return a Roadway generated from a DXF file

    Layers with names such as seg001 will contain LWPOLYLINEs.
    Each LWPOLYLINE corresponds to a lane centerline, which together
    are all neighbored.
"""
function read_dxf(io::IO, ::Type{Roadway};
    dist_threshold_lane_connect::Float64 = 0.25, # [m]
    desired_distance_between_curve_samples::Float64 = 1.0; # [m]
    )

    lines = readlines(io)

    i = findfirst(line -> contains(line, "ENTITIES"), lines)
    i != 0 || error("ENTITIES section not found")

    ###################################################
    # Pull pts for each lane
    lane_pts_dict = Dict{LaneTag, Vector{VecE2}}()

    i = findnext(lines, "LWPOLYLINE\n", i)
    while i != 0
        i = findnext(lines, "  8\n", i)
        if i != 0 # segment identifier found in LWPOLYLINE

            if ismatch(r"(?<=seg)(\d*)", lines[i+1])
                segid = parse(Int, match(r"(?<=seg)(\d*)", lines[i+1]).match)

                    i = findnext(lines, "AcDbPolyline\n", i)
                i != 0 || error("AcDbPolyline not found in LWPOLYLINE!")
                i = findnext(lines, " 90\n", i)
                i != 0 || error("Number of vertices not found in AcDbPolyline!")

                N = parse(Int, lines[i+1])
                N > 0 || error("Empty line segment!")

                pts = Array{VecE2}(N)

                i = findnext(lines, " 10\n", i)
                i != 0 || error("First point not found in AcDbPolyline!")

                for j in 1 : N
                    x = parse(Float64, lines[i+1])
                    y = parse(Float64, lines[i+3])
                    i += 4
                    pts[j] = VecE2(x,y)
                end

                laneid = 1
                for tag in keys(lane_pts_dict)
                    if tag.segment == segid
                        laneid += 1
                    end
                end
                lane_pts_dict[LaneTag(segid, laneid)] = pts
            end

            i = findnext(lines, "LWPOLYLINE\n", i)
        end
    end

    ###################################################
    # Shift pts to connect to previous / next pts
    lane_next_dict = Dict{LaneTag, Tuple{VecE2, LaneTag}}()
    lane_prev_dict = Dict{LaneTag, Tuple{VecE2, LaneTag}}()

    for (tag, pts) in lane_pts_dict
        # see if can connect to next
        best_tag = NULL_LANETAG
        best_ind = -1
        best_sq_dist = dist_threshold_lane_connect
        for (tag2, pts2) in lane_pts_dict
            if tag2.segment != tag.segment
                for (ind,pt) in enumerate(pts2)
                    sq_dist = abs2(pt - pts[end])
                    if sq_dist < best_sq_dist
                        best_sq_dist = sq_dist
                        best_ind = ind
                        best_tag = tag2
                    end
                end
            end
        end
        if best_tag != NULL_LANETAG
            # remove our last pt and set next to pt to their pt
            pop!(pts)
            lane_next_dict[tag] = (lane_pts_dict[best_tag][best_ind], best_tag)
            if best_ind == 1 # set connect prev as well
                lane_prev_dict[best_tag] = (pts[end], tag)
            end
        end
    end
    for (tag, pts) in lane_pts_dict
        # see if can connect to prev
        if !haskey(lane_prev_dict, tag)
            best_tag = NULL_LANETAG
            best_ind = -1
            best_sq_dist = dist_threshold_lane_connect
            for (tag2, pts2) in lane_pts_dict
                if tag2.segment != tag.segment
                    for (ind,pt) in enumerate(pts2)
                        sq_dist = abs2(pt - pts[1])
                        if sq_dist < best_sq_dist
                            best_sq_dist = sq_dist
                            best_ind = ind
                            best_tag = tag2
                        end
                    end
                end
            end
            if best_tag != NULL_LANETAG
                # connect 'em
                shift!(pts)
                lane_prev_dict[tag] = (lane_pts_dict[best_tag][best_ind], best_tag)
            end
        end
    end

    ###################################################
    # Build the roadway
    retval = Roadway()
    for (tag, pts) in lane_pts_dict
        if !has_segment(retval, tag.segment)
            push!(retval.segments, RoadSegment(tag.segment))
        end
    end

    lane_new_dict = Dict{LaneTag, LaneTag}() # old -> new tag
    for seg in retval.segments

        # pull lanetags for this seg
        lanetags = LaneTag[]
        for tag in keys(lane_pts_dict)
            if tag.segment == seg.id
                push!(lanetags, tag)
            end
        end

        # sort the lanes such that the rightmost lane is lane 1
        # do this by taking the first lane,
        # then project each lane's midpoint to the perpendicular at the midpoint

        @assert(!isempty(lanetags))
        proj_positions = Array{Float64}(length(lanetags))

        first_lane_pts = lane_pts_dict[lanetags[1]]
        n = length(first_lane_pts)
        lo = first_lane_pts[div(n,2)]
        hi = first_lane_pts[div(n,2)+1]
        midpt_orig = (lo + hi)/2
        dir = polar(1.0, atan2(hi - lo) + π/2) # direction perpendicular (left) of lane

        for (i,tag) in enumerate(lanetags)
            pts = lane_pts_dict[tag]
            n = length(pts)
            midpt = (pts[div(n,2)] + pts[div(n,2)+1])/2
            proj_positions[i] = proj(midpt - midpt_orig, dir, Float64)
        end

        for (i,j) in enumerate(sortperm(proj_positions))

            tag = lanetags[j]
            boundary_left = i == length(proj_positions) ? LaneBoundary(:solid, :white) : LaneBoundary(:broken, :white)
            boundary_right = i == 1 ? LaneBoundary(:solid, :white) : LaneBoundary(:broken, :white)

            pts = lane_pts_dict[tag]
            pt_matrix = Array{Float64}(2, length(pts))
            for (k,P) in enumerate(pts)
                pt_matrix[1,k] = P.x
                pt_matrix[2,k] = P.y
            end

            println("fitting curve ", length(pts), "  "); tic()
            curve = _fit_curve(pt_matrix, desired_distance_between_curve_samples)
            toc()

            tag_new = LaneTag(seg.id, length(seg.lanes)+1)
            lane = Lane(tag_new, curve,
                        boundary_left = boundary_left,
                        boundary_right = boundary_right)
            push!(seg.lanes, lane)
            lane_new_dict[tag] = tag_new
        end
    end

    ###################################################
    # Connect the lanes
    for (tag_old, tup) in lane_next_dict
        next_pt, next_tag_old = tup
        lane = retval[lane_new_dict[tag_old]]
        next_tag_new = lane_new_dict[next_tag_old]
        dest = retval[next_tag_new]
        roadproj = proj(VecSE2(next_pt, 0.0), dest, retval)

        # println("connecting $(lane.tag) to $(dest.tag)")

        cindS = curveindex_end(lane.curve)
        cindD = roadproj.curveproj.ind

        if cindD == CURVEINDEX_START # a standard connection
            connect!(lane, dest)
            # remove any similar connection from lane_prev_dict
            if haskey(lane_prev_dict, next_tag_old) && lane_prev_dict[next_tag_old][2] == tag_old
                delete!(lane_prev_dict, next_tag_old)
            end
        else
            # otherwise connect as before
            unshift!(lane.exits,  LaneConnection(true,  cindS, RoadIndex(cindD, dest.tag)))
            push!(dest.entrances, LaneConnection(false, cindD, RoadIndex(cindS, lane.tag)))
        end
    end
    for (tag_old, tup) in lane_prev_dict
        prev_pt, prev_tag_old = tup
        lane = retval[lane_new_dict[tag_old]]
        prev_tag_new = lane_new_dict[prev_tag_old]
        prev = retval[prev_tag_new]
        roadproj = proj(VecSE2(prev_pt, 0.0), prev, retval)

        # println("connecting $(lane.tag) from $(prev.tag)")

        cindS = roadproj.curveproj.ind
        cindD = CURVEINDEX_START

        if cindS == curveindex_end(prev) # a standard connection
            @assert(!has_prev(prev))
            connect!(prev, lane)
        else
            # a standard connection
            push!(prev.exits,  LaneConnection(true,  cindS, RoadIndex(cindD, lane.tag)))
            unshift!(lane.entrances, LaneConnection(false, cindD, RoadIndex(cindS, prev.tag)))
        end
    end

    retval
end

function write_dxf(io::IO, roadway::NGSIMRoadway)

    lines = open(readlines, Pkg.dir("NGSIM", "data", "template.dxf"))
    i = findfirst(line -> contains(line, "ENTITIES"), lines)
    i != 0 || error("ENTITIES section not found")

    # write out header
    for j in 1 : i
        print(io, lines[j])
    end

    # write out the lanes
    for (handle_int, lane) in enumerate(roadway.centerlines)
        write_lwpolyline(io, lane, handle_int)
    end

    # write out tail
    for j in i+1 : length(lines)
        print(io, lines[j])
    end
end
function write_roadways_to_dxf()

    roadway_input_80 = RoadwayInputParams(Pkg.dir("NGSIM", "data", "boundaries80.txt"),
                                          Pkg.dir("NGSIM", "data", "centerlines80.txt"))
    roadway_input_101 = RoadwayInputParams(Pkg.dir("NGSIM", "data", "boundaries101.txt"),
                                           Pkg.dir("NGSIM", "data", "centerlines101.txt"))

    ngsimroadway_80 = read_roadway(roadway_input_80)
    ngsimroadway_101 = read_roadway(roadway_input_101)

    open(io->write_dxf(io, ROADWAY_80), Pkg.dir("NGSIM", "data", "ngsim_80.dxf"), "w")
    open(io->write_dxf(io, ROADWAY_101), Pkg.dir("NGSIM", "data", "ngsim_101.dxf"), "w")
end
function write_roadways_from_dxf()

    roadway_80 = open(io->read_dxf(io, Roadway, dist_threshold_lane_connect=2.0), Pkg.dir("NGSIM", "data", "ngsim_80.dxf"), "r")
    roadway_101 = open(io->read_dxf(io, Roadway, dist_threshold_lane_connect=2.0), Pkg.dir("NGSIM", "data", "ngsim_101.dxf"), "r")

    # also converts to meters
    convert_curves_feet_to_meters!(roadway_80)
    convert_curves_feet_to_meters!(roadway_101)

    open(io->write(io, roadway_80), Pkg.dir("NGSIM", "data", "ngsim_80.txt"), "w")
    open(io->write(io, roadway_101), Pkg.dir("NGSIM", "data", "ngsim_101.txt"), "w")
end

const ROADWAY_80 = open(io->read(io, MIME"text/plain"(), Roadway), Pkg.dir("NGSIM", "data", "ngsim_80.txt"), "r")
const ROADWAY_101 = open(io->read(io, MIME"text/plain"(), Roadway), Pkg.dir("NGSIM", "data", "ngsim_101.txt"), "r")


