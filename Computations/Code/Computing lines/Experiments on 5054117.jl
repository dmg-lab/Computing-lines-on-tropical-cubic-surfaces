# Experiments on 5054117
# include("schlaefliwalls.jl")
# include("To be legible/5054117.jl")

# 1. Check whether hv Mot n°8 coincides with a facet of secondary cone
MotH_hv = filter(Mot -> pm.polytope.dim(visibilityConeH(Mot))< 20, MotH)
visibilityC = visibilityConeH(MotH_hv[1])
F = pm.polytope.facet(SecCone, 6)
pm.polytope.equal_polyhedra(visibilityC,F)

# 2. Check whether hv Mot n°1 coincdes with facet of secondary cone
MotD_hv = filter(Mot -> pm.polytope.dim(visibilityConeD(Mot))< 20, MotD)
Mot = MotD_hv[2]
visibilityC = visibilityConeD(Mot)
F = pm.polytope.facet(SecCone, 13)
pm.polytope.equal_polyhedra(visibilityC,F)

# 3. Each visiblity Cone is contained in a facet of the secondary cone
for Mot in MotD_hv
    i = 0
    visCone = visibilityConeD(Mot)
    while !(pm.polytope.contains(pm.polytope.facet(SecCone, i),visCone)) && i< SecCone.N_FACETS -1
        i += 1
    end
    println(i) # If i is equal to 15 the visibility cone is not contained in a facet, if i < 15 it is contained in facet(SecCone, i)
end

# 3. Check which facets contain which visibility cones
for Mot in MotD_hv
    visCone = visibilityConeD(Mot)
    println("Motif: ", Mot[1], " is contained in facets")
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:SecCone.N_FACETS -1)
    println(hyp) # prints all facets of the secondary cone containing the visibility cone
end

# 4. Check that visibility cones coincide with faces
for Mot in MotD_hv 
    visCone = visibilityConeD(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:15)
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
    for i in 1:length(V) 
        if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
    end
    V = map(i -> i-1, filter(i -> i!=0, V))
    f = pm.polytope.face(SecCone, V)
    println("Motif ", Mot[1], " coincides with a face of the visibility cone ", pm.polytope.equal_polyhedra(f, visCone))
end

# 5. Compute Hyperplane Arrangement
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotB 
    SW = SchlaefliWall(visibilityConeB(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1)end end
end
for Mot in MotC
    SW = SchlaefliWall(visibilityConeC(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1)end end
end
for Mot in MotD 
	if !(Mot in MotD_hv) 
        SW = SchlaefliWall(visibilityConeD(Mot)) 
        if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1)end end
    end
end
for Mot in MotE 
    SW = SchlaefliWall(visibilityConeH(Mot))
    if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1)end end
end
for Mot in MotH
    visCone = visibilityConeH(Mot)
	if pm.polytope.dim(visCone) == 20 
        SW = SchlaefliWall(visCone)
        if SW != [] for W in SW SWs = cat(SWs, transpose(W), dims=1)end end
    end
end
for Mot in MotJ 
    SW = SchlaefliWall(visibilityConeJ(Mot))
    if SW != [] for W in SW cat(SWs, transpose(W), dims=1)end end
end

SWs = unique(SWs, dims = 1)

# $SWs = [[0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -2, 1],
# [0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, -2, 0, 0, 0, 2, 0, 1, -1, 0],
# [0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 2, 0, -1, -1, 1],
# [0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, -1, -2, 1],
# [0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 1, 0],
# [0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, -2, 0, 0, 0, 1, 0, 1, 1, -1],
# [0, -1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 1, 0, 0, -1, 0],
# [0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 0, 1, 0, -1, 1, 0],
# [0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, -2, 0, 1, 1, -1],
# [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 1, 2, -1],
# [0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, -1, 0, 1, -1, 0],
# [0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, -2, 0, -1, 1, 0],
# [0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, -1, 0, -1, -1, 1],
# [0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 2, -1]];


HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION

# Compute number of maximal cones, f-vector
nmc = CD.N_MAXIMAL_CONES # 36
CD.F_VECTOR
# CD.N_CONES

# 6. Compute minimal number of lines
# f_normals = Matrix{Int}(CD.FACET_NORMALS)
# # import result in OSCAR
# f_normals = [0 0 0 0 1 0 0 -1 0 0 0 0 0 -1 0 0 0 0 2 -1;
# 0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 -1;
# 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -2 0 0 1;
# 1 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 0 1 0 0 0;
# 0 0 0 1 0 0 -2 0 1 0 0 0 0 0 0 0 0 0 0 0;
# 0 0 0 1 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 0 0;
# 0 0 1 -1 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0;
# 0 0 0 0 0 0 1 0 -2 1 0 0 0 0 0 0 0 0 0 0;
# 0 1 0 0 0 0 0 0 0 -1 0 -2 0 0 0 1 0 1 1 -1;
# 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 -1 1 0 0 0 0;
# 0 0 -1 0 0 0 0 0 0 1 0 1 0 0 1 -2 0 0 0 0;
# 0 -1 0 0 0 0 0 1 0 0 0 1 0 0 0 -1 0 0 0 0;
# 0 0 1 0 0 0 0 0 0 -1 0 -1 0 0 0 1 0 -1 1 0;
# 0 1 0 0 0 0 0 -1 0 1 0 -1 0 0 0 -1 0 0 1 0;
# 0 0 0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0 0 -1 1;
# 0 0 0 0 0 1 0 0 0 -1 0 -1 0 0 0 1 0 0 0 0;
# 0 0 0 0 0 0 0 -1 0 0 0 0 0 1 0 1 0 0 -1 0;
# 0 1 1 0 0 -3 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
# 0 1 0 0 0 0 0 0 0 0 0 -2 0 0 0 -1 0 1 2 -1;
# 0 1 0 0 0 0 0 0 0 -1 0 -2 0 0 0 2 0 1 -1 0;
# 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 -1 1 0;
# 0 0 0 0 1 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 0;
# 0 0 -1 0 0 0 0 0 0 1 0 1 0 0 0 -2 0 1 1 -1]
# mcones_facets = Matrix{Int}(CD.MAXIMAL_CONES_FACETS)
# mcones_facets = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0;
# 1 1 1 1 1 1 1 1 -1 1 1 1 1 1 1 1 1 1 1 1 0 0 0;
# -1 1 1 1 1 1 1 1 0 1 0 0 0 -1 0 1 1 1 0 0 1 1 1;
# 1 1 1 1 1 1 1 1 1 1 1 0 1 -1 1 1 1 1 0 0 0 0 0;
# 1 1 1 1 1 1 1 1 -1 1 0 1 -1 1 1 1 1 1 1 0 1 0 0;
# 1 1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 1 1 -1 1 0 0 0;
# 1 1 1 1 1 1 1 1 0 1 1 0 1 -1 1 1 1 1 -1 1 0 0 0;
# 1 1 1 1 1 1 1 1 1 1 0 0 -1 -1 1 1 1 1 0 0 1 0 -1;
# 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 -1 -1 0 0 0;
# 1 1 1 1 1 1 1 1 0 1 1 0 0 -1 1 1 1 1 -1 -1 0 0 0;
# 1 1 1 1 1 1 1 1 0 1 1 0 0 -1 1 1 1 1 1 -1 0 0 0;
# 1 1 1 1 1 1 1 1 -1 1 0 0 -1 -1 1 1 1 1 1 0 1 0 0;
# -1 1 1 1 1 1 1 1 0 1 1 1 0 1 0 1 1 1 -1 -1 0 1 0;
# -1 1 1 1 1 1 1 1 0 1 1 0 0 -1 0 1 1 1 -1 -1 0 1 0;
# -1 1 1 1 1 1 1 1 0 1 1 0 0 -1 0 1 1 1 1 -1 0 1 0;
# 1 1 1 1 1 1 1 1 -1 1 1 0 1 -1 1 1 1 1 1 1 0 0 0;
# -1 1 1 1 1 1 1 1 -1 1 1 0 1 -1 0 1 1 1 1 1 0 1 0;
# -1 1 1 1 1 1 1 1 -1 1 0 0 -1 -1 0 1 1 1 1 0 1 1 0;
# -1 1 1 1 1 1 1 1 1 1 1 0 1 -1 0 1 1 1 0 0 0 1 0;
# -1 1 1 1 1 1 1 1 0 1 0 0 -1 -1 0 1 1 1 -1 0 1 1 0;
# -1 1 1 1 1 1 1 1 0 1 1 0 1 -1 0 1 1 1 -1 1 0 1 0;
# -1 1 1 1 1 1 1 1 1 1 0 1 -1 1 0 1 1 1 0 0 1 1 -1;
# -1 1 1 1 1 1 1 1 0 1 1 1 0 1 0 1 1 1 1 -1 0 1 0;
# -1 1 1 1 1 1 1 1 0 1 0 1 -1 1 0 1 1 1 -1 0 1 1 0;
# -1 1 1 1 1 1 1 1 0 1 1 1 1 1 0 1 1 1 -1 1 0 1 0;
# 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 -1 0 0 0;
# 1 1 1 1 1 1 1 1 0 1 0 0 0 -1 1 1 1 1 0 0 1 0 1;
# 1 1 1 1 1 1 1 1 0 1 0 0 -1 -1 1 1 1 1 -1 0 1 0 0;
# 1 1 1 1 1 1 1 1 0 1 0 1 -1 1 1 1 1 1 -1 0 1 0 0;
# -1 1 1 1 1 1 1 1 -1 1 0 1 -1 1 0 1 1 1 1 0 1 1 0;
# -1 1 1 1 1 1 1 1 1 1 0 0 -1 -1 0 1 1 1 0 0 1 1 -1;
# -1 1 1 1 1 1 1 1 -1 1 1 1 1 1 0 1 1 1 1 1 0 1 0;
# -1 1 1 1 1 1 1 1 0 1 0 1 0 1 0 1 1 1 0 0 1 1 1;
# 1 1 1 1 1 1 1 1 0 1 0 1 0 1 1 1 1 1 0 0 1 0 1;
# 1 1 1 1 1 1 1 1 1 1 0 1 -1 1 1 1 1 1 0 0 1 0 -1;
# -1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 0 0 1 0]

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, HA)
write("SchlaefliFan5054117.json", Polymake.call_function(:common, :encode_json, serialized))

max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
r = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(CD.LINEALITY_SPACE)

for i in 1:nmc
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = r[raysid, :]
    c = Polymake.polytope.Cone(INPUT_RAYS = rays_mcone, LINEALITY_SPACE = lin_space)

# for i in 1:nrows(mcones_facets)
#     cone_facets = Vector{Int}(mcones_facets[i,:])
#     ineq_pos = f_normals[filter(j -> cone_facets[j] > 0,1:ncols(mcones_facets)), :]
#     ineq = vcat(ineq_pos,-f_normals[filter(j -> cone_facets[j] < 0, 1:ncols(mcones_facets)),:])
#     cone = pm.polytope.Cone(INEQUALITIES=ineq)

    count = 0
    for Mot in MotA
        visCone = visibilityConeA(Mot)
        check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
        if length(check) == nrows(rays_mcone) 
            count += 1 
            if !pm.polytope.equal_polyhedra(visCone,SecCone) println(Mot[1]) end 
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotB
        visCone = visibilityConeB(Mot)
        check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
        if length(check) == nrows(rays_mcone) 
            count += 1 
            if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotC
        visCone = visibilityConeC(Mot)
        check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
        if length(check) == nrows(rays_mcone) 
            count += 1 
            if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotD
        if !(Mot in MotD_hv)
            visCone = visibilityConeD(Mot)
            check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
            if length(check) == nrows(rays_mcone) 
                count += 1 
                if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
            end
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotE
        visCone = visibilityConeE(Mot)
        check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
        if length(check) == nrows(rays_mcone) 
            count += 1 
            if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotH
        if !(Mot in MotH_hv)
            visCone = visibilityConeH(Mot)
            check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
            if length(check) == nrows(rays_mcone) 
                count += 1 
                if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
            end
        end
        # if pm.polytope.contains(visCone, c) count += 1 end
    end   
    familycount = 0 
    # for Mot in MotJ
    #     visCone = visibilityConeJ(Mot)
    #     check = filter(j ->  pm.polytope.contains(visCone, rays_mcone[j,:]), 1:nrows(rays_mcone))
    #     if length(check) == nrows(rays_mcone) 
    #         count += 1 
    #         if !pm.polytope.equal_polyhedra(visCone, SecCone) println(Mot[1]) end 
    #     end
    #     # if pm.polytope.contains(visCone, c) familycount += 1 end
    # end
    println(count, ", ", familycount)
end

weights_5054117 = [12,-131,-131,-95,-67,-129,-108,-116,-115,-119,-9,-131,-92,-76,-117,-119,-24,-83,-82,-36]

for i in 1:nmc
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = r[raysid, :]
    c = Polymake.polytope.Cone(INPUT_RAYS = rays_mcone, LINEALITY_SPACE = lin_space)
    if pm.polytope.contains(c,weight_5) return i end
end
# 15
raysid = filter(j -> max_cones[15, j] != 0, 1:ncols(max_cones))
rays_mcone = r[raysid, :]
c = Polymake.polytope.Cone(INPUT_RAYS = rays_mcone, LINEALITY_SPACE = lin_space)

Mot_in_15 = []
for Mot in MotA
    visCone = visibilityConeA(Mot)
    if pm.polytope.contains(visCone, c) push!(Mot_in_15, [:A, Mot]) end
end
for Mot in MotB
    visCone = visibilityConeB(Mot)
    if pm.polytope.contains(visCone, c) push!(Mot_in_15, [:B, Mot]) end
end
for Mot in MotD
    visCone = visibilityConeD(Mot)
    if pm.polytope.contains(visCone, c) push!(Mot_in_15, [:D, Mot]) end
end
for Mot in MotH
    visCone = visibilityConeH(Mot)
    if pm.polytope.contains(visCone, c) push!(Mot_in_15, [:H, Mot]) end
end

coefficients = [c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19]
vertices_hampe = []
for Mot in Mot_in_15
    if Mot[1] == :A
        V = verticesA(Mot[2])
        V_hampe = map(j -> map(i -> evaluate(V[j][i],coefficients, weights_5054117), 1:3),1:2)
        push!(vertices_hampe, V_hampe)
    end
    if Mot[1] == :B
        V = verticesB(Mot[2])
        V_hampe = map(j -> map(i -> evaluate(V[j][i],coefficients, weights_5054117), 1:3),1:2)
        push!(vertices_hampe, V_hampe)
    end
    if Mot[1] == :D
        V = verticesD(Mot[2])
        V_hampe = map(j -> map(i -> evaluate(V[j][i],coefficients, weights_5054117), 1:3),1:2)
        push!(vertices_hampe, V_hampe)
    end
    if Mot[1] == :H
        V = verticesH(Mot[2])
        V_hampe = map(j -> map(i -> evaluate(V[j][i],coefficients, weights_5054117), 1:3),1:2)
        push!(vertices_hampe, V_hampe)
    end
end

for i in 1:nmc
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = r[raysid, :]
    c = Polymake.polytope.Cone(INPUT_RAYS = rays_mcone, LINEALITY_SPACE = lin_space)
    # v = [1,-178,-60,121,-88,-91,67,-94,21,-20,-25,-148,38,-80,-33,-72,-33,-36,-64,-31]
    v = PointVector(check)
    if pm.polytope.contains_in_interior(c, v)
        println(i)
    end
end