using Oscar
include("schlaefliwalls.jl")

function howmanylines(n::Int, vec::PointVector)
    count = 0
    visDict = Dict{Vector{Matrix{Int64}}, Polymake.BigObjectAllocated}()
    for Mot in MotA
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeA(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end
    for Mot in MotB
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeB(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end

    for Mot in MotC
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeC(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end

    for Mot in MotD
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeD(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end

    for Mot in MotE
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeE(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end

    for Mot in MotH
        if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeH(Mot)) end
        visCone = visDict[Mot]
        if pm.polytope.contains(visCone, vec) count += 1 end
    end

    fcount = 0
    if n == 5054117
       	for Mot in MotJ
            if !haskey(visDict, Mot) push!(visDict, Mot=>visibilityConeJ(Mot)) end
            visCone = visDict[Mot]
            if pm.polytope.contains(visCone, vec) fcount += 1 end
        end 
    end
    return(count, fcount)
end

println("Which triangulation?")
n = readline()
n = parse(Int64, n)

# Import motifs and seccondary cone
include("info triangulation/$n.jl")

println("What are the coefficients of f defining the surface? (Write as x0 x1 x2...)")
to_check = readline()
vec = map(i -> parse(Int, i), split(to_check))

(nlines, nfamilies) = howmanylines(n, PointVector(vec))
println("The corresponding surface has this many visible motifs of motifs 3A, 3B, 3C, 3D, 3E, 3H, i.e. lines:")
println(nlines)

if n == 5054117
    println("Attention two motifs of 3D have to be subtracted!")
    println("Further it has ", nfamilies+2, " families of lines.")
end


