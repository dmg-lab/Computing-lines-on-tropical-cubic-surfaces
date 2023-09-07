function verticesE_old(Mot)
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
    right = hcat(Dir[Mot[2][4]+1,:], Dir[Mot[2][1]+1,:])
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

function visibilityConeE_old(Mot)
    (V1, V2) = verticesE_old(Mot)
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