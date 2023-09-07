# ID: 37

# Secondary Cone
SecEq = [1 0 0 0 -2 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
0 -1 1 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0;
0 -1 0 0 1 1 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 -2 0 1 0 0 -1 1 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 -2 0 1 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 -1 0 0 1 0 0 1 0 0 -1 0 0 0 0 0;
0 0 -1 0 0 1 0 0 0 0 0 0 1 0 -1 0 0 0 0 0;
0 0 1 0 0 -1 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 1 0 -1 0 0 0 -1 0 1 0 0 0 0 0;
0 0 1 0 0 0 0 0 0 0 0 -1 -1 0 0 0 0 1 0 0;
0 0 0 0 0 1 0 -1 0 0 0 -1 0 1 0 0 0 0 0 0;
0 0 0 0 0 0 0 -1 0 1 1 0 0 -1 0 0 0 0 0 0;
0 0 0 0 0 0 0 1 0 -1 0 0 0 -1 0 1 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -1 0 1 1 -1 0 0 0 0;
0 0 0 0 0 0 0 0 0 1 0 0 0 -1 0 -1 1 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 -1 1 0 0 0 1 -1 0 0;
0 0 0 0 0 0 0 0 0 0 0 1 0 0 -2 1 -1 1 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 -1 1 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -2 1;]

SecCone = pm.polytope.Cone(INEQUALITIES=SecEq)
FSecCone = Matrix{Int64}(SecCone.FACETS)*Vars

# Motifs
MotA = [[[16 15 11 14 5 9],[3 0 2 1]],[[2 8 11 14 9 15],[1 0 2 3]],[[16 15 11 14 2 5],[3 0 2 1]],[[10 7 11 5 9 14],[3 1 2 0]],[[2 8 11 14 15 16],[1 0 2 3]]]
MotB = []
MotC = []
MotD = [[[8 14 2 5 1 7 10],[2 0 1 3]],[[9 14 11 5 1 7 10],[2 0 1 3]],[[9 15 13 11 10 1 5],[3 0 2 1]],[[14 15 16 11 10 1 5],[3 0 2 1]],[[14 17 16 11 10 1 5],[3 0 2 1]],[[7 10 1 5 2 8 14],[2 3 1 0]],[[7 9 11 5 2 8 14],[2 3 1 0]],[[7 10 11 5 2 8 14],[2 3 1 0]],[[9 15 11 14 12 2 8],[2 3 0 1]],[[15 16 11 14 12 2 8],[2 3 0 1]],[[15 16 17 14 12 2 8],[2 3 0 1]],[[7 9 5 11 2 12 14],[1 3 2 0]],[[7 10 5 11 2 12 14],[1 3 2 0]],[[9 15 13 11 10 5 7],[3 0 2 1]],[[14 15 16 11 10 5 7],[3 0 2 1]],[[14 17 16 11 10 5 7],[3 0 2 1]],[[1 5 10 11 16 14 15],[3 1 2 0]],[[5 7 10 11 16 14 15],[3 1 2 0]],[[7 9 13 11 16 14 15],[3 1 2 0]],[[1 5 10 11 16 14 17],[3 1 2 0]],[[5 7 10 11 16 14 17],[3 1 2 0]],[[7 9 13 11 16 14 17],[3 1 2 0]],[[2 5 11 14 17 15 16],[2 1 0 3]],[[5 9 11 14 17 15 16],[2 1 0 3]],[[2 8 12 14 17 15 16],[2 1 0 3]]]
MotE = [[[12 11 14 9 15 5 9],[0 2 3 1]],[[17 11 14 9 15 5 9],[0 2 3 1]],[[12 11 14 15 16 5 9],[0 2 3 1]],[[17 11 14 15 16 5 9],[0 2 3 1]],[[12 11 14 9 15 2 5],[0 2 3 1]],[[17 11 14 9 15 2 5],[0 2 3 1]],[[12 11 14 15 16 2 5],[0 2 3 1]],[[17 11 14 15 16 2 5],[0 2 3 1]],[[1 5 11 9 14 7 10],[1 2 0 3]],[[2 5 11 9 14 7 10],[1 2 0 3]],[[1 5 11 9 14 7 9],[1 2 0 3]],[[2 5 11 9 14 7 9],[1 2 0 3]],[[10 11 13 9 15 7 9],[2 3 0 1]],[[16 11 13 9 15 7 9],[2 3 0 1]]]
MotH = []
MotJ = []

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