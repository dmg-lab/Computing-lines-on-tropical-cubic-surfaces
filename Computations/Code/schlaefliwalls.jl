# Ring 
PC = [0 0 0; 0 0 1; 0 0 2; 0 0 3;0 1 0;0 1 1;0 1 2;0 2 0;0 2 1;0 3 0;1 0 0;1 0 1;1 0 2;1 1 0;1 1 1;1 2 0;2 0 0;2 0 1;2 1 0;3 0 0]
R, (x, y, z, a, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19) = QQ["x", "y", "z", "a", "c0","c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13","c14","c15","c16","c17","c18","c19"]
S33 = matrix_space(R, 3, 3)
S23 = matrix_space(R, 2, 3)
S200 = matrix_space(R, 0, 20)

# Variables c0,...,c19
Vars = gens(R)[5:24]

# Tropical Monomials
FF = [c0; c1+z; c2+2*z; c3+3*z; c4+y; c5+y+z; c6+y+2*z; c7+2*y; c8+2*y+z; c9+3*y; c10+x; c11+x+z; c12+x+2*z; c13+x+y; c14+x+y+z; c15+x+2*y; c16+2*x; c17+2*x+z; c18+2*x+y; c19+3*x]

# Exits 
Dir = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1]

# Secondary cone
function matrix_from_vectors(vectors)
	M = vectors[1]
	for i in 2:length(vectors)
		M = hcat(M, vectors[i])
	end
	return matrix(QQ, M)
end
# Sec = [c0-2*c10 + c16; c1 + c2 - 3*c5 + c9; c4-2*c13+c18; c10-2*c16+c19; -c4 + c13+c16-c19; c5-c9-c11+c15; -c1+c7+c11-c15;  -c7 + c13 + c15 - c18; c7 - c13 - c18 + c19; c3 - 2*c6 + c8; c2-c3-c11+c12; c6-2*c8+c9; c8-c9-c14+c15; -c2 + c9 + c11 + c14 - 2*c15; c3 - 2*c12 + c17; c14 - c15 - c17 + c18]
# SecEq = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Sec)))
# SecCone = cone_from_inequalities(SecEq)
# FSecCone = facets(SecCone)

function pts_in_motif(points)
	Mpts = PC[points[1]+1:points[1]+1, :]
	for i = 2:length(points)
		Mpts = vcat(Mpts, PC[points[i]+1:points[i]+1, :])
	end
	return Mpts
end	

function coeff_in_motif(points)
	Mcoef = Vars[points[1]+1]
	for i = 2:length(points)
		Mcoef = vcat(Mcoef, Vars[points[i]+1])
	end
	return Mcoef
end

function matrix_from_vectors(vectors)
	M = vectors[1]
	for i in 2:length(vectors)
		M = hcat(M, vectors[i])
	end
	return matrix(QQ, M)
end

# Polymake
const pm = Polymake

# Method to compute vertices of lines dual to occurrences of motif A. 
function verticesA(Mot) #Mot is a Vector 2 elements, 1. points of Motif, 2. combinatorial type
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])
	
	# Vertex CDEF
	M1 = matrix(R, vcat(Mpts[3:3,:] - Mpts[4:4,:], vcat(Mpts[4:4,:] - Mpts[5:5,:], Mpts[5:5,:] - Mpts[6:6,:])))
	N1 = [Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]; Mcoef[6] - Mcoef[5]]
	VertCDEF = inv(M1)N1

	# Vertex ACD
	Ray = VertCDEF - S33(a)*(Dir[Mot[2][4] + 1,:])
	M2 = matrix(R, vcat(Mpts[1:1,:] - Mpts[3:3,:], Mpts[3:3,:] - Mpts[4:4,:]))
	N2 = [Mcoef[3] - Mcoef[1]; Mcoef[4] - Mcoef[3]]
	Sys2 = M2*Ray - N2
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	VertACD = map(f -> evaluate(f, [a], [l]), Ray)

	# Vertex ABD
	Edge = VertACD + S33(a)*(Dir[Mot[2][1]+1, :] + Dir[Mot[2][2]+1, :])
	M3 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], Mpts[2:2,:] - Mpts[4:4,:]))
	N3 = [Mcoef[2] - Mcoef[1]; Mcoef[4] - Mcoef[2]];
	Sys3 = M3*Edge - N3;
	if Sys3[1] != 0
		if coeff(Sys3[1], a) == 1 l = -Sys3[1] + a else  l = Sys3[1] + a end
	else 
		if coeff(Sys3[2], a) == 1 l = -Sys3[2] + a else  l = Sys3[2] + a end
	end
	VertABD =  map(f -> evaluate(f, [a], [l]), Edge)
	return [VertABD, VertACD]
end

function visibilityConeA(Mot)
	(V1, V2) = verticesA(Mot)
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][1] + 1], Hyp2)
	Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Totineq)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])
	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function verticesB(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	M1 = matrix(R, vcat(Mpts[2:2,:] - Mpts[3:3,:], vcat(Mpts[3:3,:] - Mpts[4:4,:], Mpts[4:4,:] - Mpts[5:5,:])))
	N1 = [Mcoef[3] - Mcoef[2]; Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]]
	VertBCDE = inv(M1)*N1
	
	# EdgeL = VertBCDE - a * (Dir[Mot[2][1] + 1,:] + Dir[Mot[2][2] + 1])
	# EdgeR = VertBCDE + a * (Dir[Mot[2][1] + 1,:] + Dir[Mot[2][2] + 1])

	EdgeL = VertBCDE - S33(a) * (Dir[Mot[2][1] + 1, :] + Dir[Mot[2][2] + 1,:])
	M2 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], Mpts[2:2,:] - Mpts[3:3,:]))
	N2 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]]
	Sys2 = M2 * EdgeL - N2
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	VertABC = map(f -> evaluate(f, [a], [l]), EdgeL)

	EdgeR = VertBCDE + S33(a) * (Dir[Mot[2][1] + 1, :] + Dir[Mot[2][2] + 1, :])
	M3 = matrix(R, vcat(Mpts[4:4,:] - Mpts[5:5,:], Mpts[5:5,:] - Mpts[6:6,:]))
	N3 = [Mcoef[5] - Mcoef[4]; Mcoef[6] - Mcoef[5]]
	Sys3 = M3 * EdgeR - N3
	if Sys3[1] != 0
		if coeff(Sys3[1], a) == 1 l = -Sys3[1] + a else  l = Sys3[1] + a end
	else
		if coeff(Sys3[2], a) == 1 l = -Sys3[2] + a else  l = Sys3[2] + a end
	end
	VertDEF = map(f -> evaluate(f, [a],  [l]), EdgeR)
	return [VertABC, VertDEF]
end

function visibilityConeB(Mot)
	(V1, V2) = verticesB(Mot)
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][4] + 1], Hyp2)
	Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Totineq)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])
	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function verticesC(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	# Vertex BCDE
	M1 = matrix(R, vcat(Mpts[2:2,:] - Mpts[3:3,:], vcat(Mpts[3:3,:] - Mpts[4:4,:], Mpts[4:4,:] - Mpts[5:5,:])));
	N1 = [Mcoef[3] - Mcoef[2]; Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]]
	VertBCDE = inv(M1)*N1 

	# Vertex DEFG
	M2 = matrix(R, vcat(Mpts[4:4,:] - Mpts[5:5,:], vcat(Mpts[5:5,:] - Mpts[6:6,:], Mpts[6:6,:] - Mpts[7:7,:])));
	N2 = [Mcoef[5] - Mcoef[4]; Mcoef[6] - Mcoef[5]; Mcoef[7] - Mcoef[6]]
	VertDEFG = inv(M2)*N2 

	# Vertex DE = intrsection of VertBCDE - l*(Dir[Mot[2][1]+1,:] + Dir[Mot[2][2]+1,:]) and VertDEFG - m*Dir[Mot[2][4]+1,:], l, m > 0
	left = VertBCDE - VertDEFG
	right = hcat(Dir[Mot[2][1]+1,:] + Dir[Mot[2][2]+1,:], Dir[Mot[2][4]+1,:])
	(i, j) = (1,2)
	while det(vcat(right[i:i,:], right[j:j,:])) == 0
		if j<=2 j += 1; else i += 1 end
	end
	(l,m) = inv(matrix(R, vcat(right[i:i,:], right[j:j,:]))) * vcat(left[i], left[j])
	VertDE = VertDEFG - S33(m)*Dir[Mot[2][4]+1 ,:]

	Ray = VertBCDE + S33(a)*(Dir[Mot[2][1]+1,:] + Dir[Mot[2][2]+1,:])
	M3 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:],Mpts[2:2,:] - Mpts[3:3,:])) 
	N3 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]]
	Sys2= M3 * Ray - N3
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	VertABC = map(f -> evaluate(f, [a], [l]), Ray)
	return [VertABC, VertDE]
end

function visibilityConeC(Mot)
	(V1, V2) = verticesC(Mot)
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	
	Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][4] + 1], Hyp2)
	Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Totineq)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])

	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function verticesD(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	M1 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], vcat(Mpts[2:2,:] - Mpts[3:3,:], Mpts[3:3,:] - Mpts[4:4,:])));
	N1 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]; Mcoef[4] - Mcoef[3]]
	VertABCD = inv(M1)*N1

	Ray1 = VertABCD - S33(a)* Dir[Mot[2][2]+1, :]
	M2 = matrix(R, vcat(Mpts[3:3,:] - Mpts[4:4,:],Mpts[4:4,:] - Mpts[5:5,:]))
	N2 = [Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]];
	Sys2= M2 * Ray1 - N2;
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	VertCDE = map(f -> evaluate(f, [a], [l]), Ray1)

	M3 = matrix(R, vcat(Mpts[4:4,:] - Mpts[5:5,:], vcat(Mpts[5:5,:] - Mpts[6:6,:], Mpts[6:6,:] - Mpts[7:7,:])));
	N3 = [Mcoef[5] - Mcoef[4]; Mcoef[6] - Mcoef[5]; Mcoef[7] - Mcoef[6]]
	VertDEFG = inv(M3)*N3;
	  
	Ray2 = VertDEFG - S33(a)*Dir[Mot[2][4] + 1, :]
	right = hcat(Dir[Mot[2][1] + 1, :] + Dir[Mot[2][2] + 1, :], Dir[Mot[2][4] + 1, :])
	left = VertDEFG - VertCDE;
	(i, j) = (1,2)
	while det(vcat(right[i:i,:], right[j:j,:])) == 0
		if j <= 2 j += 1 else i += 1 end
	end
	(l,m) = inv(matrix(R, right[i:j,:])) * vcat(left[i], left[j]);
	Ver2 = map(f -> evaluate(f, [a], [m]), Ray2)
	return (VertCDE, Ver2)
end

function visibilityConeD(Mot)
	# Vertices of line
	(V1, V2) = verticesD(Mot)
	# The minimum has to be attained in cells corresponding to motif
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	Ineq1 = map(f -> f - Hyp1[Mot[1][3] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][4] + 1], Hyp2)
	Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Totineq)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])
	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function verticesE(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	# Vertex BCDE
	M1 = matrix(R, vcat(Mpts[2:2,:] - Mpts[3:3,:], vcat(Mpts[3:3,:] - Mpts[4:4,:], Mpts[4:4, :] - Mpts[5:5,:])));
	N1 = [Mcoef[3] - Mcoef[2]; Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]]
	VertBCDE = inv(M1)*N1

	# Vertex BCFG
	M2 = matrix(R, vcat(Mpts[2:2,:] - Mpts[3:3,:], vcat(Mpts[3:3,:] - Mpts[6:6,:], Mpts[6:6, :] - Mpts[7:7,:])))
	N2 = [Mcoef[3] - Mcoef[2]; Mcoef[6] - Mcoef[3]; Mcoef[7] - Mcoef[6]]
	VertBCFG = inv(M2)*N2

	# Vertex BC
	left = VertBCDE - VertBCFG
	right = hcat(Dir[Mot[2][4]+1,:], -Dir[Mot[2][3]+1,:])
	(i, j) = (1,2)
	while det(vcat(right[i:i,:], right[j:j,:])) == 0
		if j <= 2 j += 1; else i += 1 end
	end
	(l,k) = inv(matrix(R, vcat(right[i:i,:], right[j:j,:]))) * vcat(left[i], left[j])
	VertBC = VertBCDE - S33(l)*Dir[Mot[2][4]+1,:]

	Ray = VertBC + S33(a)*(Dir[Mot[2][1]+1, :] + Dir[Mot[2][2]+1,:])
	M3 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], Mpts[2:2,:] - Mpts[3:3,:]))
	N3 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]]
	Sys1 = M3*Ray - N3
	if Sys1[1] != 0
		if coeff(Sys1[1], a) == 1 l = -Sys1[1] + a else  l = Sys1[1] + a end
	else
		if coeff(Sys1[2], a) == 1 l = -Sys1[2] + a else  l = Sys1[2] + a end
	end
	VertABC = map(f -> evaluate(f, [a], [l]), Ray)
	return [VertABC, VertBC]
end

function visibilityConeE(Mot)
	(V1, V2) = verticesE(Mot)
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	
	Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][3] + 1], Hyp2)
	Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Totineq)))
	c = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])

	inter = pm.polytope.intersection(c, SecCone)
	return inter
end


function verticesH(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	M1 = matrix(R, vcat(Mpts[1:1,:]-Mpts[2:2,:], vcat(Mpts[2:2,:]-Mpts[3:3,:],Mpts[3:3,:]-Mpts[4:4,:])))
	N1 = [Mcoef[2]-Mcoef[1];Mcoef[3]-Mcoef[2];Mcoef[4]-Mcoef[3]]
	VertABCD = inv(M1)*N1

	Edge = VertABCD + S33(a)*(Dir[Mot[2][1]+1,:] + Dir[Mot[2][2]+1,:])
	M2 = matrix(R, vcat(Mpts[3:3,:] - Mpts[4:4,:],Mpts[4:4,:] - Mpts[5:5,:]))
	N2 = [Mcoef[4] - Mcoef[3]; Mcoef[5] - Mcoef[4]];
	Sys2 = M2*Edge - N2;
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	VertCDE = map(f -> evaluate(f, [a], [l]), Edge)
	return [VertABCD, VertCDE]
end	

function visibilityConeH(Mot)
	(V1, V2) = verticesH(Mot)
	# Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Hyp2 = map(f-> evaluate(f, [x, y, z], [V2[1], V2[2], V2[3]]), FF);
	# Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	Ineq2 = map(f -> f - Hyp2[Mot[1][5] + 1], Hyp2)
	# Totineq = vcat(Ineq1, Ineq2)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Ineq2)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])
	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function verticesJ(Mot)
	Mpts = pts_in_motif(Mot[1])
	Mcoef = coeff_in_motif(Mot[1])

	M1 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], vcat(Mpts[2:2,:] - Mpts[3:3,:], Mpts[3:3,:] - Mpts[4:4,:])))
	N1 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]; Mcoef[4] - Mcoef[3]]
	VertABCD = inv(M1)*N1

	# M2 = matrix(R, vcat(Mpts[1:1,:] - Mpts[2:2,:], vcat(Mpts[2:2,:] - Mpts[3:3,:], Mpts[3:3,:] - Mpts[5:5,:])))
	# N2 = [Mcoef[2] - Mcoef[1]; Mcoef[3] - Mcoef[2]; Mcoef[5] - Mcoef[3]]
	# VertABCE = inv(M2)*N2

	Ray1 = VertABCD - S33(a)*Dir[Mot[2][2]+1,:]
	# Ray2 = VertABCE - S33(a)*Dir[Mot[2][1]+1,:]

	M3 = matrix(R, vcat(Mpts[1:1,:]-Mpts[4:4,:],Mpts[4:4,:]-Mpts[5:5,:]))
	N3 = [Mcoef[4]-Mcoef[1]; Mcoef[5]-Mcoef[4]]
	Sys2 = M3 * Ray1 - N3
	if Sys2[1] != 0
		if coeff(Sys2[1], a) == 1 l = -Sys2[1] + a else  l = Sys2[1] + a end
	else
		if coeff(Sys2[2], a) == 1 l = -Sys2[2] + a else  l = Sys2[2] + a end
	end
	return map(f -> evaluate(f, [a], [l]), Ray1)
end

function visibilityConeJ(Mot)
	V1 = verticesJ(Mot)
	Hyp1 = map(f-> evaluate(f, [x, y, z], [V1[1], V1[2], V1[3]]), FF);
	Ineq1 = map(f -> f - Hyp1[Mot[1][1] + 1], Hyp1)
	IM = transpose(matrix_from_vectors(map(f -> map(v -> coeff(f, v), Vars), Ineq1)))
	cone = pm.polytope.Cone(INEQUALITIES=IM[filter(i->!iszero(IM[i,:]), 1:nrows(IM)), :])
	inter = pm.polytope.intersection(cone, SecCone)
	return inter
end

function SchlaefliWall(visibilityC)
	FVisCone = Matrix{Int64}(visibilityC.FACETS)*Vars
	SchlaefliW = filter(f -> !(f in FSecCone), FVisCone)
	if pm.polytope.dim(visibilityC) < 20
		Equations = Matrix{Rational{Int64}}(visibilityC.LINEAR_SPAN)*Vars
		println("Equations: ", Equations)
	end
	SchlaefliWEq = map(f -> map(v -> Int(coeff(f, v)), Vars), SchlaefliW)
	return SchlaefliWEq
end

# function SchlaefliWall(visibilityC)
# 	FVisCone = Matrix{Int64}(visibilityC.FACETS)*Vars
# 	if pm.polytope.dim(visibilityC) < 20
# 		hyp = filter(i -> pm.polytope.contains(pm.polytope.facet(SecCone,i),visibilityC), 0:Int(SecCone.N_FACETS)-1)
# 		V = Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[1]+1, :])
# 		for i in 2:length(hyp) V+= Vector{Int}(SecCone.RAYS_IN_FACETS[hyp[i]+1,:]) end
# 		for i in 1:length(V) 
# 			if V[i]==length(hyp) V[i] = i else V[i] = 0 end 
# 		end
# 		V = map(i -> i-1, filter(i -> i!=0, V))
# 		f = pm.polytope.face(SecCone, V)
# 		if !pm.polytope.equal_polyhedra(f, visibilityC)
# 			Equations = Matrix{Rational{Int64}}(visibilityC.LINEAR_SPAN)*Vars
# 			FSecE = Matrix{Int64}(pm.polytope.intersection(SecCone, visibilityC.LINEAR_SPAN).FACETS)*Vars
# 			SchlaefliW = filter(f -> !(f in FSecE), FVisCone)
# 			return (SchlaefliW, Equations)
# 		else
# 			return []
# 		end
# 	end
# 	SchlaefliW = filter(f -> !(f in FSecCone), FVisCone)
# 	SchlaefliWEq = map(f -> map(v -> Int(coeff(f, v)), Vars), SchlaefliW)
# 	return SchlaefliWEq
# end
