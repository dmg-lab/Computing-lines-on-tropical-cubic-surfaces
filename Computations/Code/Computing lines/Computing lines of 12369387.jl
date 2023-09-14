# Experiments on triangulation #12359387 (honeycomb triangulation)
using Oscar
include("../schlaefliwalls.jl")
include("../info triangulation/12369387.jl")

# 1. Compute Schläfli fan
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotB
    SW = SchlaefliWall(visibilityConeB(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotD 
	SW = SchlaefliWall(visibilityConeD(Mot)) 
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotE
    SW = SchlaefliWall(visibilityConeE(Mot))
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end

SWs = unique(SWs, dims = 1)
HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
CD = HA.CHAMBER_DECOMPOSITION
nmc = CD.N_MAXIMAL_CONES 

# serialize and save Schläfli fan
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, CD)
write("SchlaefliFan12369387.json", Polymake.call_function(:common, :encode_json, serialized))

# 3. Compute minimal number of lines
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
    println(count + 1) # +1 for globally visible motifs 3F
    push!(counts, count)
end

println("Number(s) of lines:")
println(counts)
