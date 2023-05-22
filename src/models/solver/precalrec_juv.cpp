#include "calpop.h"

///Forward function for:
///precalrec for larval and juvenile life stages. 
///This routine precomputes diagonal coefficient for calrec_adre 

void CCalpop::precalrec_juv_comp(const PMap& map, dmatrix& bm, const dmatrix& mortality)
{
	//bm = b;
	for (int j = b.rowmin(); j <= b.rowmax(); j++){
		for (int i = b[j].indexmin(); i <= b[j].indexmax(); i++){
			bm[j][i] += mortality[i][j];
		}
	}
} 
