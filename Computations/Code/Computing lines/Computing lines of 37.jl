# Experiments on triangulation #37
using Oscar
include("../schlaefliwalls.jl")
include("../info triangulation/37.jl")

# 0. Preliminaries
MotE_hv = filter(Mot -> pm.polytope.dim(visibilityConeE(Mot)) != 20, MotE) # == MotE

# 1. Check which hyperplanes of secondary cone contain Motif 
hypDict = Dict{Vector{Matrix{Int}}, Vector{Int}}()

for Mot in MotE_hv
    visCone = visibilityConeE(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:18)
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
    for i in 1:length(V) 
        if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
    end
    V = map(i -> i-1, filter(i -> i!=0, V))
    f = pm.polytope.face(SecCone, V)
    println("The visibility cone of motif ", Mot[1], " is contained in a face of the secondary cone: ", pm.polytope.contains(f, visCone))
    println("The visibility cone of motif ", Mot[1], " coincides with a face of the secondary cone: ", pm.polytope.equal_polyhedra(f, visCone))
    println("It lies in hyperplanes", hyp)
    hypDict[Mot] = hyp
end
println(" ")

# 2. Check if two following visibility cones form a face of the secondary cone
i = 1
lin_space = SecCone.LINEALITY_SPACE
println("The visibility cones of hardly visible motifs i and i+1 coincide, where i is an even number:")
while i < length(MotE)
    visCone_1 = visibilityConeE(MotE[i])
    visCone_2 = visibilityConeE(MotE[i+1])
    union_cone = pm.polytope.Cone(INPUT_RAYS=vcat(visCone_1.RAYS, visCone_2.RAYS), INPUT_LINEALITY=lin_space)
    hyp = hypDict[MotE[i]]
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for j in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[j]+1,:]) end
    for j in 1:length(V) 
        if V[j]==length(hyp) V[j] = j else V[j] = 0 end 
    end
    V = map(j -> j-1, filter(j -> j!=0, V))
    f = pm.polytope.face(SecCone, V)
    println(pm.polytope.equal_polyhedra(f, union_cone))
    global i += 2
end
println(" ")

# 3. Compute Hyperplane Arrangement
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotD 
	SW = SchlaefliWall(visibilityConeD(Mot)) 
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims= 1) end end
end
SWs = unique(SWs, dims = 1)
HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION
nmc = CD.N_MAXIMAL_CONES # 199

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, CD)
write("SchlaefliFan37.json", Polymake.call_function(:common, :encode_json, serialized))

# 4. Compute minimal number of lines
max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
r = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(CD.LINEALITY_SPACE)

visDict = Dict{Vector{Matrix{Int64}}, Polymake.BigObjectAllocated}()
counts = Set{Int}()

for i in 1:nmc
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = r[raysid, :]
    c = pm.polytope.Cone(INPUT_RAYS=rays_mcone, INPUT_LINEALITY=lin_space)
    count = 0
    for Mot in MotA
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeA(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end
    
    for Mot in MotB
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeB(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end

    for Mot in MotC
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeC(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end

    for Mot in MotD
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeD(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end

    for Mot in MotE
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeE(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end

    for Mot in MotH
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeH(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) count += 1 end
    end 
    push!(counts, count + 1) # +1 for globally visible motifs 3F
end

println("Number(s) of lines:")
println(counts)

