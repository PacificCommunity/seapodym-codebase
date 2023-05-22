#include "calpop.h"

///Forward main function called in simulation mode only for: 
///precalrec for larval and juvenile life stages. 
///See precalrec_juv.cpp

void CCalpop::Precalrec_juv(const PMap& map, CMatrices& mat, dvar_matrix& mortality, const int t_count)
{
	dmatrix M_c    = value(mortality);
	dmatrix bm_c   = value(dvarsBM);
	dmatrix xbet_c = value(Xbet);

	dvar_matrix W(0,maxn-1,0,maxn-1);
	W.initialize();

	bm_c = b;

	precalrec_juv_comp(map, bm_c, M_c);

	dvarsBM = nograd_assign(bm_c);

	bm = bm_c;
	
	xbet_comp(map, xbet_c, a, bm_c, c, 2*iterationNumber);
	
	Xbet = nograd_assign(xbet_c);
	
}
