# ID: 5054117

# Secondary Cone
SecEq = [1 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 0 1 0 0 0;
0 1 1 0 0 -3 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 1 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 0;
0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 -2 0 0 1;
0 0 0 0 -1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 -1;
0 0 0 0 0 1 0 0 0 -1 0 -1 0 0 0 1 0 0 0 0;
0 -1 0 0 0 0 0 1 0 0 0 1 0 0 0 -1 0 0 0 0;
0 0 0 0 0 0 0 -1 0 0 0 0 0 1 0 1 0 0 -1 0;
0 0 0 0 0 0 0 1 0 0 0 0 0 -1 0 0 0 0 -1 1;
0 0 0 1 0 0 -2 0 1 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 -1 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0;
0 0 0 0 0 0 1 0 -2 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 1 -1 0 0 0 0 -1 1 0 0 0 0;
0 0 -1 0 0 0 0 0 0 1 0 1 0 0 1 -2 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 -1 1 0]

M = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1]

SecCone = pm.polytope.Cone(INEQUALITIES=SecEq)
FSecCone = Matrix{Int64}(SecCone.FACETS)*Vars

# Motifs
MotA = [[[18 17 15 11 2 9],[0 2 3 1]],[[18 19 15 11 2 9],[3 2 0 1]],[[18 19 15 11 2 9],[0 2 3 1]],[[18 17 15 11 1 9],[0 2 3 1]],[[18 19 15 11 1 9],[3 2 0 1]],[[18 19 15 11 1 9],[0 2 3 1]]]
MotB = [[[9 15 7 1 18 19],[0 1 2 3]],[[17 18 11 1 15 7],[0 2 1 3]],[[17 18 11 1 15 9],[0 2 1 3]],[[19 18 11 1 15 7],[0 2 1 3]],[[19 18 11 1 15 9],[0 2 1 3]]]
MotC = []
MotD = [[[3 14 2 11 1 9 15],[1 0 2 3]],[[9 15 2 11 1 9 15],[1 0 2 3]],[[14 15 2 11 1 9 15],[1 0 2 3]],[[3 14 2 11 1 15 18],[1 0 2 3]],[[9 15 2 11 1 15 18],[1 0 2 3]],[[14 15 2 11 1 15 18],[1 0 2 3]],[[3 14 2 11 1 18 19],[1 0 2 3]],[[9 15 2 11 1 18 19],[1 0 2 3]],[[14 15 2 11 1 18 19],[1 0 2 3]],[[9 15 1 11 2 3 14],[1 3 2 0]],[[15 18 1 11 2 3 14],[1 3 2 0]],[[18 19 1 11 2 3 14],[1 3 2 0]],[[9 15 1 11 2 9 15],[1 3 2 0]],[[15 18 1 11 2 9 15],[1 3 2 0]],[[18 19 1 11 2 9 15],[1 3 2 0]],[[9 15 1 11 2 14 15],[1 3 2 0]],[[15 18 1 11 2 14 15],[1 3 2 0]],[[18 19 1 11 2 14 15],[1 3 2 0]],[[2 3 14 11 17 15 18],[0 1 2 3]],[[1 9 15 11 17 15 18],[0 1 2 3]],[[2 9 15 11 17 15 18],[0 1 2 3]],[[2 3 14 11 17 18 19],[0 1 2 3]],[[1 9 15 11 17 18 19],[0 1 2 3]],[[2 9 15 11 17 18 19],[0 1 2 3]]]
MotE = []
MotH = [[[9 15 2 14 3],[1 3 2 0]],[[18 19 1 13 4],[0 2 1 3]],[[7 15 1 18 19],[0 1 2 3]],[[18 19 1 13 7],[0 2 1 3]],[[9 11 1 15 7],[0 2 1 3]],[[11 18 1 15 7],[0 2 1 3]],[[11 18 1 15 9],[0 2 1 3]]]
MotJ = [[[11 9 15 1 2],[0 3 1 2]]]

for Mot in MotA
	println("MotifA: ", Mot[1])
	visibilityC = visibilityConeA(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotB 
	println("MotifB: ", Mot[1])
	visibilityC = visibilityConeB(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotC
	println("MotifC: ", Mot[1])
	visibilityC = visibilityConeC(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotD 
	println("MotifD: ", Mot[1])
	visibilityC = visibilityConeD(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotE 
	println("MotifE: ", Mot[1])
	visibilityC = visibilityConeE(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotH
	println("MotifH: ", Mot[1])
	visibilityC = visibilityConeH(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end
for Mot in MotJ 
	println("MotifJ: ", Mot[1])
	visibilityC = visibilityConeJ(Mot)
	SWs = SchlaefliWall(visibilityC)
	for SW in SWs println(transpose(Vars)*SW) end
end

