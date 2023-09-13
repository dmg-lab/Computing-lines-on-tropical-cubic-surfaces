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