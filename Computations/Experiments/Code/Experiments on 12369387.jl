# 1. Compute Schläfli fan
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotB
    SW = SchlaefliWall(visibilityConeB(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotD 
	SW = SchlaefliWall(visibilityConeD(Mot)) 
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotE
    SW = SchlaefliWall(visibilityConeE(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
end

SWs = unique(SWs, dims = 1)

HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION

nmc = CD.N_MAXIMAL_CONES 

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, CD)
write("SchlaefliFan12369387.json", Polymake.call_function(:common, :encode_json, serialized))

# 3. Compute minimal number of lines
# f_normals = Matrix{Int}(CD.FACET_NORMALS)
# mcones_facets = Matrix{Int}(CD.MAXIMAL_CONES_FACETS)
# f_vector = CD.F_VECTOR
max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
rays = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(CD.LINEALITY_SPACE)

count = 0
for i in 1:nmc
    # cone_facets = mcones_facets[i,:]
    # ineq_pos = f_normals[filter(i -> cone_facets[i] > 0,1:ncols(mcones_facets)),:]
    # ineq = vcat(ineq_pos,-f_normals[filter(i -> cone_facets[i] < 0,1:ncols(mcones_facets)),:])
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = rays[raysid, :]
    c = pm.polytope.Cone(INPUT_RAYS=rays_mcone, INPUT_LINEALITY=lin_space)
    count = 0
    for Mot in MotA
        visCone = visibilityConeA(Mot)
        if pm.polytope.contains(visCone, c) count += 1 end
    end
    for Mot in MotB
        visCone = visibilityConeB(Mot)
        if pm.polytope.contains(visCone, c) count += 1 end
    end
    for Mot in MotD
        visCone = visibilityConeD(Mot)
        if pm.polytope.contains(visCone, c) count += 1 end
    end
    for Mot in MotE
        visCone = visibilityConeE(Mot)
        if pm.polytope.contains(visCone, c) count += 1 end
    end
    if count != 19 println("Attention: ", i) end
end
