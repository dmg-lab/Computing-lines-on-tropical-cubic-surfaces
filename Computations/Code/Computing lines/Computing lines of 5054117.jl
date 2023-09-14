# Experiments on 5054117
using Oscar
include("../schlaefliwalls.jl")
include("../info triangulation/5054117.jl")

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

# 4. Compute Hyperplane Arrangement
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotB 
    SW = SchlaefliWall(visibilityConeB(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotC
    SW = SchlaefliWall(visibilityConeC(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotD 
	if !(Mot in MotD_hv) 
        SW = SchlaefliWall(visibilityConeD(Mot)) 
        if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
    end
end
for Mot in MotE 
    SW = SchlaefliWall(visibilityConeH(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1 )end end
end
for Mot in MotH
    visCone = visibilityConeH(Mot)
	if !(Mot in MotH_hv)
        SW = SchlaefliWall(visCone)
        if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
    end
end
for Mot in MotJ 
    SW = SchlaefliWall(visibilityConeJ(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end

SWs = unique(SWs, dims = 1)
HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION
nmc = CD.N_MAXIMAL_CONES # 36
# CD.F_VECTOR # Takes very long!

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, HA)
write("SchlaefliFan5054117.json", Polymake.call_function(:common, :encode_json, serialized))

# 5. Compute minimal number of lines
max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
r = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(CD.LINEALITY_SPACE)

visDict = Dict{Vector{Matrix{Int64}}, Polymake.BigObjectAllocated}()
counts = Dict{Symbol, Set{Int}}(:isolated => Set{Int}(), :family => Set{Int}())

for i in 1:nmc
    raysid = filter(j -> max_cones[i, j] != 0, 1:ncols(max_cones))
    rays_mcone = r[raysid, :]
    c = Polymake.polytope.Cone(INPUT_RAYS = rays_mcone, LINEALITY_SPACE = lin_space)
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

    familycount = 0
    for Mot in MotJ
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeJ(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, c) familycount += 1 end
    end
    println(count + 4, ", ", familycount+2) # +4 for 6 globally visible motifs 3F, 3G and -2 for two motifs 3D, +2 for globally visible motifs 3I
    push!(counts, :isolated => push!(counts[:isolated], count+4))
    push!(counts, :family => push!(counts[:family], familycount+2))
end

println("Number of isolated lines:")
println(counts[:isolated])
println("Number of family of lines:")
println(counts[:family])