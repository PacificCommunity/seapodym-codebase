#include "calpop.h"

void CCalpop::InitCalPop(CParam& param, const PMap& map) 	
{
	a.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	b.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	bm.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	c.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	d.allocate(map.imin, map.imax, map.jinf, map.jsup);
	e.allocate(map.imin, map.imax, map.jinf, map.jsup);
	f.allocate(map.imin, map.imax, map.jinf, map.jsup);
	xbet.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	ybet.allocate(map.imin, map.imax, map.jinf, map.jsup);
	a.initialize();
	b.initialize();
	bm.initialize();
	c.initialize();
	d.initialize();
	e.initialize();
	f.initialize();
	xbet.initialize();
	ybet.initialize();

	dvarsA.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	dvarsB.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	dvarsBM.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	dvarsC.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	dvarsD.allocate(map.imin, map.imax, map.jinf, map.jsup);
	dvarsE.allocate(map.imin, map.imax, map.jinf, map.jsup);
	dvarsF.allocate(map.imin, map.imax, map.jinf, map.jsup);
	Xbet.allocate(map.jmin, map.jmax, map.iinf, map.isup);
	Ybet.allocate(map.imin, map.imax, map.jinf, map.jsup);
	dvarsA.initialize();
	dvarsB.initialize();
	dvarsBM.initialize();
	dvarsC.initialize();
	dvarsD.initialize();
	dvarsE.initialize();
	dvarsF.initialize();
	Xbet.initialize();
	Ybet.initialize();

 	const int nti = param.get_nbi();
	const int ntj = param.get_nbj();
	uuint.allocate(0, nti, 0, ntj);
	uuint.initialize();

	maxn = Utilities::MyMax( nti, ntj );//param.get_maxn();
	iterationNumber = param.iterationNumber; 	
	deltax = param.deltaX;
	deltay = param.deltaY;
	sigma_fcte = param.sigma_fcte;


	//The maximal speed, which is supported by ADI, i.e. does not generate negative 
	//densities due to approximation errors with finite differences in case of strong 
	//gradients, high MSS (mostly for largest fish) and nearly zero densities. Vinf 
	//will be used to constrain resulting velocities as a function of local habitat 
	//gradients. Vinf is set to allow <=nmax grid cells moves per one iteration. The
	//value nmax is set empirically. Thus, from optimizations with SKJ tagging data
	//it is necessary to set nmax=2.
	//Examples of run configurations with 
	//nmax=3: 
	//for 1deg,1month,30N Vinf = 5400 nmi/mo; for 1/4deg,7days,28N Vinf = 5400 nmi/mo
	//for 1/12deg,1day,12N Vinf = 5400 nmi/mo
	//nmax=2: 
	//for 1deg,1month,45N Vinf = 5400 nmi/mo; for 1/4deg,7days,42N Vinf = 5400 nmi/mo
	//for 1/12deg,1day,18N Vinf = 5400 nmi/mo
	const int nmax = 1;
	Vinf = nmax * deltax * iterationNumber; 

	//Exploited biomass
	int nb_ne_fishery = sum(param.mask_fishery_sp_no_effort[0]);//ATTN: violation of multi-species
	cout << nb_ne_fishery << endl;
	dvarsSNsum.allocate(0,nb_ne_fishery-1);
	for (int fne=0; fne<nb_ne_fishery; fne++){
		dvarsSNsum(fne).allocate(map.imin, map.imax, map.jinf, map.jsup);
		dvarsSNsum(fne).initialize();
	}

	//Selectivity
	const int nb_fishery = param.get_nbfishery();
	const int nb_species = param.get_nbspecies();
	Selectivity.allocate(0,nb_species-1);
	for (int sp=0; sp<nb_species; sp++){
		const int a0 = param.sp_a0_adult[sp];
		const int nb_ages = param.sp_nb_cohorts[sp];
		Selectivity(sp).allocate(0,nb_fishery-1,a0,nb_ages-1);
		Selectivity(sp).initialize();	
	}
}
