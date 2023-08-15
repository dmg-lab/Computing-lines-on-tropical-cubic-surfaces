# ID: 842554

# Secondary Cone
SecEq = [1 -1 0 0 0 0 0 0 0 0 -1 1 0 0 0 0 0 0 0 0;
0 -1 1 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 -1 -1 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 -2 0 1 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 -1 0 0 0 0 -1 0 0 1 0 0 0 0 0;
0 0 -1 0 0 1 1 0 -1 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 -1 0 0 -1 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 -1 0 0 1 0 0 0 0 0 1 0 -1 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 0 0;
0 0 0 0 0 1 0 -1 0 0 0 -1 0 1 0 0 0 0 0 0;
0 0 0 0 1 0 0 -2 0 1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 -1 0 0 1 0 0 1 0 0 -1 0 0 0 0 0 0;
0 0 0 0 1 0 0 0 0 0 -1 0 0 -1 0 0 1 0 0 0;
0 0 0 0 0 -1 0 1 1 -1 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 1 0 -1 0 0 0 -1 0 1 0 0 0 0;
0 0 0 0 0 0 0 0 1 -1 0 0 0 0 -1 1 0 0 0 0;
0 0 0 0 0 0 0 0 0 1 0 0 0 -1 0 -1 1 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 -2 0 1;
0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 -1 -1 0 1 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 -2 1;
0 0 0 0 0 0 0 0 0 0 0 -1 0 0 1 0 1 0 -1 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 1 1 -1]

SecCone = pm.polytope.Cone(INEQUALITIES=SecEq)
FSecCone = Matrix{Int64}(SecCone.FACETS)*Vars

# Motifs
MotA = [[[4 10 5 11 8 9],[3 2 1 0]],[[18 16 14 11 8 9],[3 2 0 1]],[[18 19 14 11 8 9],[3 2 0 1]],[[18 16 14 11 3 6],[3 2 0 1]],[[18 19 14 11 3 6],[3 2 0 1]],[[4 1 13 11 9 15],[1 2 3 0]],[[6 2 14 11 9 15],[1 2 0 3]],[[6 3 14 11 9 15],[1 2 0 3]],[[6 2 14 11 15 18],[1 2 0 3]],[[6 3 14 11 15 18],[1 2 0 3]],[[6 2 14 11 18 19],[1 2 0 3]],[[6 3 14 11 18 19],[1 2 0 3]],[[18 16 14 11 6 8],[3 2 0 1]],[[18 19 14 11 6 8],[3 2 0 1]]]
MotB = []
MotC = []
MotD = [[[9 15 13 11 10 1 4],[3 0 2 1]],[[15 18 16 11 10 1 4],[3 0 2 1]],[[18 19 16 11 10 1 4],[3 0 2 1]],[[3 6 2 11 1 4 10],[1 0 2 3]],[[6 8 2 11 1 4 10],[1 0 2 3]],[[8 9 5 11 1 4 10],[1 0 2 3]],[[4 10 1 11 2 3 6],[1 3 2 0]],[[4 7 5 11 2 3 6],[1 3 2 0]],[[7 9 5 11 2 3 6],[1 3 2 0]],[[4 10 1 11 2 6 8],[1 3 2 0]],[[4 7 5 11 2 6 8],[1 3 2 0]],[[7 9 5 11 2 6 8],[1 3 2 0]],[[1 4 10 11 16 15 18],[3 1 2 0]],[[4 7 13 11 16 15 18],[3 1 2 0]],[[7 9 13 11 16 15 18],[3 1 2 0]],[[1 4 10 11 16 18 19],[3 1 2 0]],[[4 7 13 11 16 18 19],[3 1 2 0]],[[7 9 13 11 16 18 19],[3 1 2 0]]]
MotE = [[[3 11 14 9 15 8 9],[0 2 3 1]],[[12 11 14 9 15 8 9],[0 2 3 1]],[[17 11 14 9 15 8 9],[0 2 3 1]],[[19 11 14 9 15 8 9],[0 2 3 1]],[[3 11 14 15 18 8 9],[0 2 3 1]],[[12 11 14 15 18 8 9],[0 2 3 1]],[[17 11 14 15 18 8 9],[0 2 3 1]],[[19 11 14 15 18 8 9],[0 2 3 1]],[[3 11 14 18 19 8 9],[0 2 3 1]],[[12 11 14 18 19 8 9],[0 2 3 1]],[[17 11 14 18 19 8 9],[0 2 3 1]],[[19 11 14 18 19 8 9],[0 2 3 1]],[[3 11 14 9 15 3 6],[0 2 3 1]],[[12 11 14 9 15 3 6],[0 2 3 1]],[[17 11 14 9 15 3 6],[0 2 3 1]],[[19 11 14 9 15 3 6],[0 2 3 1]],[[3 11 14 15 18 3 6],[0 2 3 1]],[[12 11 14 15 18 3 6],[0 2 3 1]],[[17 11 14 15 18 3 6],[0 2 3 1]],[[19 11 14 15 18 3 6],[0 2 3 1]],[[3 11 14 18 19 3 6],[0 2 3 1]],[[12 11 14 18 19 3 6],[0 2 3 1]],[[17 11 14 18 19 3 6],[0 2 3 1]],[[19 11 14 18 19 3 6],[0 2 3 1]],[[1 5 11 8 9 7 9],[1 2 0 3]],[[2 5 11 8 9 7 9],[1 2 0 3]],[[10 11 13 9 15 7 9],[2 3 0 1]],[[16 11 13 9 15 7 9],[2 3 0 1]],[[1 5 11 8 9 4 7],[1 2 0 3]],[[2 5 11 8 9 4 7],[1 2 0 3]],[[10 11 13 9 15 4 7],[2 3 0 1]],[[16 11 13 9 15 4 7],[2 3 0 1]],[[3 11 14 9 15 6 8],[0 2 3 1]],[[12 11 14 9 15 6 8],[0 2 3 1]],[[17 11 14 9 15 6 8],[0 2 3 1]],[[19 11 14 9 15 6 8],[0 2 3 1]],[[3 11 14 15 18 6 8],[0 2 3 1]],[[12 11 14 15 18 6 8],[0 2 3 1]],[[17 11 14 15 18 6 8],[0 2 3 1]],[[19 11 14 15 18 6 8],[0 2 3 1]],[[3 11 14 18 19 6 8],[0 2 3 1]],[[12 11 14 18 19 6 8],[0 2 3 1]],[[17 11 14 18 19 6 8],[0 2 3 1]],[[19 11 14 18 19 6 8],[0 2 3 1]]]
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

