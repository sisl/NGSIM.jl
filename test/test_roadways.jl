p = lerp(CurvePt(VecSE2(0.0,0.0,0.0), 0.0), CurvePt(VecSE2(1.0,2.0,3.0), 4.0), 0.25)
@test isapprox(convert(VecE2, p.pos), VecE2(0.25, 0.5))
@test isapprox(p.pos.θ, 0.75)
@test isapprox(p.s, 1.0)

@test isapprox(_mod2pi2(0.0), 0.0)
@test isapprox(_mod2pi2(0.5), 0.5)
@test isapprox(_mod2pi2(2pi + 1), 1.0)
@test isapprox(_mod2pi2(1 - 2pi), 1.0)

centerline = [CurvePt(VecSE2(0.0,0.0,0.0), 0.0),
			  CurvePt(VecSE2(1.0,0.0,0.0), 1.0),
			  CurvePt(VecSE2(3.0,0.0,0.0), 3.0)]
@test _binary_search_curve_dist2(centerline, VecE2(0.0,0.0)) == 1
@test _binary_search_curve_dist2(centerline, VecE2(1.0,0.0)) == 2
@test _binary_search_curve_dist2(centerline, VecE2(2.1,0.0)) == 3
@test _binary_search_curve_dist2(centerline, VecE2(0.49,0.0)) == 1
@test _binary_search_curve_dist2(centerline, VecE2(1.9,-100.0)) == 2
@test _binary_search_curve_dist2(centerline, VecSE2(1.9,-100.0,0.0)) == 2
@test _binary_search_curve_dist2(centerline, VecSE2(-1.0,0.0,0.0)) == 1

@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(1.0,0.0), VecE2(1.0,0.0)), 1.0)
@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(1.0,0.0), VecE2(0.5,0.0)), 0.5)
@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(1.0,0.0), VecE2(0.5,0.5)), 0.5)
@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(1.0,1.0), VecE2(0.5,0.5)), 0.5)
@test isapprox(_proj_rel(VecE2(1.0,0.0), VecE2(2.0,0.0), VecE2(1.5,0.5)), 0.5)
@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(-1.0,0.0), VecE2(1.0,0.0)), 0.0)
@test isapprox(_proj_rel(VecE2(0.0,0.0), VecE2(-1.0,0.0), VecE2(-0.75,0.0)), 0.75)

@test _get_ind_lo_and_hi(centerline, 1.0) == (1,2)
@test _get_ind_lo_and_hi(centerline, 2.5) == (2,3)
@test _get_ind_lo_and_hi(centerline, 3.0) == (2,3)

@test isapprox(get_extind(centerline, 0.0), 1.0)
@test isapprox(get_extind(centerline, 1.0), 2.0)
@test isapprox(get_extind(centerline, 1.5), 2.25)
@test isapprox(get_extind(centerline, 3.5), 3.0)
@test isapprox(get_extind(centerline,-0.5), 1.0)

@test isapprox(move_extind_along(1.0, centerline, 0.25), 1.25)
@test isapprox(move_extind_along(1.0, centerline, 1.25), 2.125)
@test isapprox(move_extind_along(1.5, centerline, 1.5), 2.25)

res = project_to_lane(VecSE2(0.0,0.0,0.0), centerline)
@test isapprox(res.extind, 1.0)
@test isapprox(res.t, 0.0)
@test isapprox(res.ϕ, 0.0)
res = project_to_lane(VecSE2(0.25,0.5,0.1), centerline)
@test isapprox(res.extind, 1.25)
@test isapprox(res.t, 0.5)
@test isapprox(res.ϕ, 0.1)
res = project_to_lane(VecSE2(0.25,-0.5,-0.1), centerline)
@test isapprox(res.extind, 1.25)
@test isapprox(res.t, -0.5)
@test isapprox(res.ϕ, -0.1)
res = project_to_lane(VecSE2(1.5,0.5,-0.1), centerline)
@test isapprox(res.extind, 2.25)
@test isapprox(res.t, 0.5)
@test isapprox(res.ϕ, -0.1)

centerline2 = [CurvePt(VecSE2(-1.0,-1.0,π/4), 0.0),
			   CurvePt(VecSE2( 1.0, 1.0,π/4), hypot(2.0,2.0)),
			   CurvePt(VecSE2( 3.0, 3.0,π/4), hypot(4.0,4.0))]

res = project_to_lane(VecSE2(0.0,0.0,0.0), centerline2)
@test isapprox(res.extind, 1.5)
@test isapprox(res.t, 0.0)
@test isapprox(res.ϕ, -π/4)

centerlines = Array(Vector{CurvePt}, 2)
centerlines[1] = centerline
centerlines[2] = centerline2
testroad = Roadway(:test, Vector{VecE2}[], centerlines)

@test isapprox(project_posG_to_frenet(VecSE2(0.0,0.0,0.0), testroad), Frenet(1,1.0,0.0,0.0,0.0))
@test isapprox(project_posG_to_frenet(VecSE2(3.0,3.0,0.0), testroad), Frenet(2,3.0,hypot(4.0,4.0),0.0,-π/4))
@test isapprox(project_posG_to_frenet(VecSE2(1.5,-1.0,-0.1), testroad), Frenet(1,2.25,1.5,-1.0,-0.1))