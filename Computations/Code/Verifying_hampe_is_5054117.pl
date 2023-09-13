$F = toTropicalPolynomial("min(12+3*x0,-131+2*x0+x1, -67+2*x0+x2,-9+2*x0+x3,-131+x0+2*x1,-129+x0+x1+x2, -131+x0+x1+x3,-116+x0+2*x2,-76+x0+x2+x3,-24+x0+2*x3,-95+3*x1, -108+2*x1+x2,-92+2*x1+x3,-115+x1+2*x2,-117+x1+x2+x3, -83+x1+2*x3,-119+3*x2,-119+2*x2+x3,-82+x2+2*x3,-36+3*x3)");
application "tropical";
$V = new Hypersurface<Min>(POLYNOMIAL=>$F);
$mcF = new Set<Set>(rows($V->DUAL_SUBDIVISION->MAXIMAL_CELLS));

application "fan";
$X = retrieve_by_id(5054117);
$mcX = new Set<Set>(rows($X->MAXIMAL_CELLS));

application "group";
$Xpts = new Matrix($X->POINTS->minor(All, ~[0]));
$Vpts = new Matrix($V->DUAL_SUBDIVISION->POINTS->minor(All,[2,3,4]));
$perm = find_permutation(rows($Vpts), rows($Xpts));
$Mperm = new Matrix(permutation_matrix($perm));
print $Mperm; 
# copy $Mperm into julia and multiply with $F->coefficients_as_vector; 
# the result is:
$weights = [12,-9,-24,-36,-67,-76,-82,-116,-119,-119,-131,-131,-83,-129,-117,-115,-131,-92,-108,-95];
$dualF = new fan::SubdivisionOfPoints(POINTS=>$X->POINTS, WEIGHTS=>$weights);
$mcF = new Set<Set>(rows($dualF->MAXIMAL_CELLS));
$o = orbit<on_elements>($X->ACTION, $mcX);
print $o->contains($mcF);
# true
$sym = [0,10,16,19,4,13,18,7,15,9,1,11,17,5,14,8,2,12,6,3];
$Msym = new Matrix(permutation_matrix($sym));
# copy $Msym into julia and multiply with $weights; 
# the result is the representative of the orbit of $weights under S_4 in the Sec:
$weights_5054117 = [12,-131,-131,-95,-67,-129,-108,-116,-115,-119,-9,-131,-92,-76,-117,-119,-24,-83,-82,-36];
# Check whether $weights_5054117 lies in a hyperplane
