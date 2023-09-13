# Experiments on 842554
using Oscar
include("schlaefliwalls.jl")
include("To be legible/842554.jl")

MotE_hv = filter(Mot -> pm.polytope.dim(visibilityConeE(Mot))< 20, MotE)

# 1. Check if Motifs are contained in facet of secondary cone
hypDict = Dict{Vector{Matrix{Int}}, Vector{Int}}()
for Mot in MotE_hv 
    visCone = visibilityConeE(Mot)
    hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visCone), 0:21)
    V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
    for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
    for i in 1:length(V) 
        if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
    end
    V = map(i -> i-1, filter(i -> i!=0, V))
    hypDict[Mot] = V
    f = pm.polytope.face(SecCone, V)
    println("Motif ", Mot[1], " coincides with a face of the visibility cone ", pm.polytope.equal_polyhedra(f, visCone))
    println("Motif ", Mot[1], " is contained in a face of the visibility cone ", pm.polytope.contains(f, visCone))
end

# 2. Check if two cones form one
i = 1
lin_space = SecCone.LINEALITY_SPACE
while i < length(MotE_hv[9:16])
    visCone_1 = visibilityConeE(MotE_hv[8+i])
    visCone_2 = visibilityConeE(MotE_hv[9+i])
    union_cone = pm.polytope.Cone(INPUT_RAYS=vcat(visCone_1.RAYS, visCone_2.RAYS), INPUT_LINEALITY=lin_space)
    V = hypDict[MotE_hv[i+8]]
    W = hypDict[MotE_hv[i+9]]
    f = pm.polytope.face(SecCone, V)
    println(pm.polytope.equal_polyhedra(f, union_cone))
    i += 2
end

# 2. Check which visibility cones coincide
for i in 1:length(MotE_hv)-1
    visCone = visibilityConeE(MotE_hv[i])
    j = 1
    while !pm.polytope.equal_polyhedra(visCone, visibilityConeE(MotE_hv[i+j])) && j < length(MotE_hv)-i
        j += 1
    end
    if pm.polytope.equal_polyhedra(visCone, visibilityConeE(MotE_hv[i+j])) println(i-1, ", ", i+j-1) end
end

# Compute Schläfli hyperplanes
SWs = Matrix{Int}(undef, 0, 20)
for Mot in MotA
    SW = SchlaefliWall(visibilityConeA(Mot))
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
nmc = CD.N_MAXIMAL_CONES
# nc = CD.N_CONES

# Compute minimal number of lines
# f_normals = Matrix{Int}(CD.FACET_NORMALS)
# mcones_facets = Matrix{Int}(CD.MAXIMAL_CONES_FACETS)

max_cones = Matrix{Int}(CD.MAXIMAL_CONES)
rays = Matrix{Rational}(CD.RAYS)
lin_space = Matrix{Rational}(SecCone.LINEALITY_SPACE)

serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, HA)
write("SchlaefliFan842554_new.json", Polymake.call_function(:common, :encode_json, serialized))

# for i in 1:nmc
#     cone_facets = mcones_facets[i,:]
    # ineq_pos = f_normals[filter(i -> cone_facets[i] > 0,1:ncols(mcones_facets)),:]
    # ineq = vcat(ineq_pos,-f_normals[filter(i -> cone_facets[i] < 0,1:ncols(mcones_facets)),:])
    # cone = pm.polytope.Cone(INEQUALITIES=ineq)
n = 842554
# println("n: ", n)
# input = "SchlaefliFan$n.json"
# mpha_result = "schlaefli_mc_$n.xz"
# println("Files: $input $mpha_result")
# arrangement = load(input)
counts = Set{Int}()
visDict = Dict()
i = 1

open(`xzcat schlaefli_mc_$n.xz`) do io
    while !eof(io)
        m = match(r"Signature: (\[.*\]) Rays: (\[\[.*\]\])", readline(io))
        if m != nothing
            sig = Meta.eval(Meta.parse(m[1]))
            R = replace(m[2], r"/" => "//")
            R = matrix(QQ, Meta.eval(Meta.parse(R)))
            c = pm.polytope.Cone(INPUT_RAYS=R, INPUT_LINEALITY=lin_space)
            count = 0
            for Mot in MotA
                if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeA(Mot)) end
                visCone = visDict[Mot]
                if pm.polytope.contains(visCone, c) count += 1 end
            end
            for Mot in MotD
                if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeD(Mot)) end
                visCone = visDict[Mot]
                if pm.polytope.contains(visCone, c) count += 1 end
            end
            for Mot in MotE
                if !(Mot in MotE_hv)
                    if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeE(Mot)) end
                    visCone = visDict[Mot]
                    if pm.polytope.contains(visCone, c) count += 1 end
                end
            end
            println(count)
            push!(counts, count)
        end
    end
end
open("counts", "w") do io
    print(io, counts)
end
println(counts)

n = 842554
n_lines = Set{Int}([14, 15, 17, 18, 19, 20, 21])
visDict = Dict()
cones = Dict{Int64, Polymake.BigObjectAllocated}()
open(`xzcat schlaefli_mc_$n.xz`) do io
    while !eof(io)
        m = match(r"Signature: (\[.*\]) Rays: (\[\[.*\]\])", readline(io))
        if m != nothing
            sig = Meta.eval(Meta.parse(m[1]))
            R = replace(m[2], r"/" => "//")
            R = matrix(QQ, Meta.eval(Meta.parse(R)))
            c = pm.polytope.Cone(INPUT_RAYS=R, INPUT_LINEALITY=lin_space)
            count = 0
            for Mot in MotA
                if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeA(Mot)) end
                visCone = visDict[Mot]
                if pm.polytope.contains(visCone, c) count += 1 end
            end
            for Mot in MotD
                if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeD(Mot)) end
                visCone = visDict[Mot]
                if pm.polytope.contains(visCone, c) count += 1 end
            end
            for Mot in MotE
                if !(Mot in MotE_hv)
                    if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeE(Mot)) end
                    visCone = visDict[Mot]
                    if pm.polytope.contains(visCone, c) count += 1 end
                end
            end
            if !haskey(cones, count) 
                push!(cones, count => c)
                pop!(n_lines, count)
            end
        end
        if is_empty(n_lines) return cones end
    end
end