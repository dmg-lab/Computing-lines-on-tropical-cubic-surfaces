MotE_hv = filter(Mot -> pm.polytope.dim(visibilityConeE(Mot)) < 20, MotE)

# 1. Check if visibility cones of hv motifs are contained in a face of secondary cone, if so check if they coincide
for Mot in MotE_hv 
    visCone = visibilityConeE(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:15)
#     println(hyp, ", ", pm.polytope.dim(visCone))
# end
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
    for i in 1:length(V) 
        if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
    end
    V = map(i -> i-1, filter(i -> i!=0, V))
    f = pm.polytope.face(SecCone, V)
    println("The visibility cone of motif ", Mot[1], " coincides with a face of the secondary cone ", pm.polytope.equal_polyhedra(f, visCone))
end

# 2. Compute Schläfli fan
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
    if !(Mot in MotE_hv)
         SW = SchlaefliWall(visibilityConeE(Mot))
        if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
    end
end

SWs = unique(SWs, dims = 1)

HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION
nmc = CD.N_MAXIMAL_CONES #26232

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, HA)
Polymake.call_function(:common, :encode_json, serialized)
write("SchlaefliFan12369387.json", Polymake.call_function(:common, :encode_json, serialized))

# 3. Compute minimal number of lines
# f_normals = Matrix{Int}(CD.FACET_NORMALS)
# mcones_facets = Matrix{Int}(CD.MAXIMAL_CONES_FACETS)
# f_vector = CD.F_VECTOR
max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
rays = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(CD.LINEALITY_SPACE)

for i in 1:nmc
    # cone_facets = mcones_facets[i,:]
    # ineq_pos = f_normals[filter(i -> cone_facets[i] > 0,1:ncols(mcones_facets)),:]
    # ineq = vcat(ineq_pos,-f_normals[filter(i -> cone_facets[i] < 0,1:ncols(mcones_facets)),:])
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = rays[raysid, :]
    cone = pm.polytope.Cone(INEQUALITIES=ineq)
    count = 0
    for Mot in MotA
        visCone = visibilityConeA(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    end
    for Mot in MotB
        visCone = visibilityConeB(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    for Mot in MotD
        visCone = visibilityConeD(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    end
    for Mot in MotE
        visCone = visibilityConeE(Mot)
        if pm.polytope.contains(visCone, cone) count += 1 end
    end
    println(count)
end
