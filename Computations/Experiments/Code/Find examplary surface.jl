r = hcat(Int.(ones(74,1)), vcat(Int.(zero(rand(1,20))), Matrix{Rational}(c.RAYS)))
r = Matrix{Rational}(r)
poly = pm.polytope.Polytope(POINTS=r)
obj = rand(Int,20)
milp = mixed_integer_linear_program(polyhedron(poly),obj)
lin_vector = Int.(Matrix{Rational}(c.LINEALITY_SPACE)[1,:])
qqzero = optimal_solution(milp)
while optimal_solution(milp) == lin_vector || optimal_solution(milp) == -lin_vector || optimal_solution(milp) == qqzero
    obj = rand(Int, 20)
    milp = mixed_integer_linear_program(polyhedron(poly),obj)
end
println(optimal_solution(milp))