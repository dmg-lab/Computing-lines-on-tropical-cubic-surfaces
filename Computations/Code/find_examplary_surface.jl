prep_r = vcat(Int.(zero(rand(21))), vcat(Matrix{Rational}(r), vcat(Matrix{Rational}(SecCone.LINEALITY_SPACE), -Matrix{Rational}(SecCone.LINEALITY_SPACE))))
r = Matrix{Rational}(hcat(ones(Int, nrows(prep_r)),  Matrix{Rational}(prep_r)))

# set n = # of triangulation
poly = pm.polytope.Polytope(POINTS=r)
poly_julia = polyhedron(-1*Matrix{Int}(c.FACETS), Int.(zero(rand(1,c.N_FACETS)))[1,:])
min_weights = Dict(12369387 => [5,1,1,5,2,0,2,2,2,5,2,0,2,0,0,1,2,2,1,5], 5054117=>[44,0,1,15,21,0,9,2,4,0,38,0,15,16,4,1,33,16,14,29], 842554 =>[10,3,1,5,3,0,2,1,0,0,6,0,6,3,2,3,7,9,8,14], 2091566 => [7,6,6,11,1,2,5,0,0,1,2,0,8,0,4,2,4,10,5,13], 37 =>[9,3,4,16,4,0,9,0,3,0,2,0,7,1,2,2,4,7,8,15])
check = min_weights[n]
b = map(bi -> trunc(Int, bi), 5*maximum(check)*rand(21))
box = polyhedron(vcat(identity_matrix(QQ, 21),-identity_matrix(QQ, 21)), vcat(b,b))
while !contains(box, check)
    b = 10*maximum(check)*rand(21)
    box = polyhedron((vcat(identity_matrix(QQ, 21),-identity_matrix(QQ, 21)), vcat(b,b)), (identity_matrix(QQ, 21)[1,:], check[21]))
end
inter = intersect(poly_julia, box)

# milp
objective = rand(-1:1,21)
milp = mixed_integer_linear_program(inter,objective)
found = false
while !found
    objective = rand(20)
    milp = mixed_integer_linear_program(inter,objective)
    if howmanylines(n, optimal_solution(milp)+check) == (21, 0) 
        println(optimal_solution(milp)+check)
        found = true
    end
end
serialized = Polymake.call_function(Symbol("Core::Serializer"), :serialize, optimal_solution(milp))
write("ex_37.json", Polymake.call_function(:common, :encode_json, serialized))
println(optimal_solution(milp))


# 842554
j = 21
include("schlaefliwalls.jl")
include("info triangulation/842554.jl")
include("compute_number_lines.jl")
c21 = load("cones_842554/n_21_lines.json")
c21_poly = polyhedron(facets(c21))
b = 500*ones(Int, 20)
box = polyhedron(vcat(identity_matrix(QQ, 20),-identity_matrix(QQ, 20)), vcat(b,b))
inter_21 = intersect(c21_poly, box)

found = false
while !found
    objective = rand(20)
    milp = mixed_integer_linear_program(inter_21, objective)
    o = optimal_solution(milp)
    if howmanylines(842554, o) == (21, 0) 
        println(o) 
        found = true
    end
end