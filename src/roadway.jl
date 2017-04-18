const FLOATING_POINT_REGEX = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
const METERS_PER_FOOT = 0.3048

type NGSIMRoadway
    name::Symbol
    boundaries::Vector{Vector{VecE2}}
    centerlines::Vector{Vector{CurvePt}}
end
immutable RoadwayInputParams
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

    retval = Array(Vector{VecE2}, n_boundaries)
    line_index = 2
    for i in 1 : n_boundaries

        @assert(lines[line_index+=1] == @sprintf("BOUNDARY %d", i))
        npts = parse(Int, lines[line_index+=1])
        line = Array(VecE2, npts)
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
        line = Array(VecE2, npts)
        for j in 1 : npts
            matches = matchall(FLOATING_POINT_REGEX, lines[line_index+=1])
            x = parse(Float64, matches[1]) * METERS_PER_FOOT # convert to meters
            y = parse(Float64, matches[2]) * METERS_PER_FOOT # convert to meters
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

function write_dxf(io::IO, roadway::NGSIMRoadway)

    lines = open(readlines, Pkg.dir("NGSIM", "data", "template.dxf"))
    i = findfirst(lines, "ENTITIES\n")
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


