# Experiments on triangulation #2091566
using Oscar
include("../schlaefliwalls.jl")
include("../info triangulation/2091566.jl")

# 0. Preliminaries
MotC_hv = filter(Mot -> pm.polytope.dim(visibilityConeC(Mot))< 20, MotC)
MotE_hv = filter(Mot -> pm.polytope.dim(visibilityConeE(Mot))< 20, MotE)

# 1. Compute faces of secondary cone containing lower dimensional visibility cones and check if they coincide
for Mot in MotC_hv 
    visCone = visibilityConeC(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:Int(SecCone.N_FACETS)-1)
    println(hyp, ", ", pm.polytope.dim(visCone))
    if hyp != []
        V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
        for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
        for i in 1:length(V) 
            if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
        end
        V = map(i -> i-1, filter(i -> i!=0, V))
        f = pm.polytope.face(SecCone, V)
        println("The visibility cone of motif ", Mot[1], " coincides with a face of the secondary cone ", pm.polytope.equal_polyhedra(f, visCone))
    end
end

for Mot in MotE_hv 
    visCone = visibilityConeE(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:Int(SecCone.N_FACETS)-1)
    println(hyp, ", ", pm.polytope.dim(visCone))
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
    for i in 1:length(V) 
        if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
    end
    V = map(i -> i-1, filter(i -> i!=0, V))
    f = pm.polytope.face(SecCone, V)
    println("The visibility cone of motif ", Mot[1], " coincides with a face of the secondary cone ", pm.polytope.contains(f, visCone))
end

# 2. Compute Schläfli hyperplanes
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
    if !(Mot in MotC_hv)
        SW = SchlaefliWall(visibilityConeC(Mot))
        if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
    end
end
for Mot in MotD 
	SW = SchlaefliWall(visibilityConeD(Mot)) 
    if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
end
for Mot in MotE
    if !(Mot in MotE_hv)
        SW = SchlaefliWall(visibilityConeE(Mot))
        if SW != [] for W in SW global SWs = cat(SWs, transpose(W), dims=1) end end
    end
end
SWs = unique(SWs, dims = 1)

HA = pm.fan.HyperplaneArrangement(HYPERPLANES=SWs, SUPPORT=SecCone)
# CD = HA.CHAMBER_DECOMPOSITION # too big to be computed sequentially, has thus been computed in mpha

# Serialze and save Schläfli fan for mpha
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, HA)
write("SchlaefliFan2091566.json", Polymake.call_function(:common, :encode_json, serialized))

# 3. Compute number of lines on a surface generic enough
lin_space = SecCone.LINEALITY_SPACE
visDict = Dict()
counts = Set{Int}()
open(`xzcat ../../Schlaefli_fans/SchlaefliFan2091566.out.xz`) do io
    while !eof(io)
        m = match(r"Signature: (\[.*\]) Rays: (\[\[.*\]\])", readline(io))
        if m != nothing
            sig = Meta.eval(Meta.parse(m[1]))
            R = replace(m[2], r"/" => "//")
            R = matrix(QQ, Meta.eval(Meta.parse(R)))
            c = pm.polytope.Cone(INPUT_RAYS=R, INPUT_LINEALITY=lin_space)
            push!(co, c)
            return 
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
            for Mot in MotD
                if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeD(Mot)) end
                visCone = visDict[Mot]
                if pm.polytope.contains(visCone, c) count += 1 end
            end
            println(count)
            push!(counts, count)
        end
    end
end

println("Number(s) of lines:")
println(counts)