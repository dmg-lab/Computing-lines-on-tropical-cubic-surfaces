MotC_hv = filter(Mot -> pm.polytope.dim(visibilityConeC(Mot))< 20, MotC)
MotE_hv = filter(Mot -> pm.polytope.dim(visibilityConeE(Mot))< 20, MotE)

SWs = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW SWs = vcat(SWs, transpose(W)) end end
end
for Mot in MotB
    SW = SchlaefliWall(visibilityConeB(Mot))
    if SW != [] for W in SW SWs = vcat(SWs, transpose(W)) end end
end
for Mot in MotC
    if !(Mot in MotC_hv)
        SW = SchlaefliWall(visibilityConeC(Mot))
        if SW != [] for W in SW SWs = vcat(SWs, transpose(W)) end end
    end
end
for Mot in MotD 
	SW = SchlaefliWall(visibilityConeD(Mot)) 
    if SW != [] for W in SW SWs = vcat(SWs, transpose(W)) end end
end
for Mot in MotE
    if !(Mot in MotE_hv)
        SW = SchlaefliWall(visibilityConeE(Mot))
        if SW != [] for W in SW SWs = vcat(SWs, transpose(W)) end end
    end
end
HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs[2:nrows(SWs),:], SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION
CD.N_MAXIMAL_CONES

# Compute minimal number of lines
f_normals = Matrix{Int}(CD.FACET_NORMALS)
mcones_facets = Matrix{Int}(CD.MAXIMAL_CONES_FACETS)

for i in 1:nrows(mcones_facets)
    cone_facets = mcones_facets[i,:]
    ineq_pos = f_normals[filter(i -> cone_facets[i] > 0,1:ncols(mcones_facets)),:]
    ineq = vcat(ineq_pos,-f_normals[filter(i -> cone_facets[i] < 0,1:ncols(mcones_facets)),:])
    cone = pm.polytope.Cone(INEQUALITIES=ineq)
    count = 0
    for Mot in MotA
        visCone = visibilityConeA(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    end
    for Mot in MotD
        visCone = visibilityConeD(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    end
    for Mot in MotE
        if !(Mot in MotE_hv)
            visCone = visibilityConeE(Mot)
            if pm.polytope.contains(visCone, cone) count += 1 end
        end
    end
    println(count)
end