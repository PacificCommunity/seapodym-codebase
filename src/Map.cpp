#include "Map.h"

////////////////////////////////////////////////////////////
// fonction lit_map de la classe PMap 
// lit le fichier source et charge les donnees dans les tableaux

int nbl;
double deltaX;
//CBord bord_layer;


void PMap::lit_map(CParam &param)
{
	const int nti=param.get_nbi();
	const int ntj=param.get_nbj();
	nbl=param.nb_layer;

// creation des tableaux


	bord_cell.allocate(0, nti - 1, 0, ntj - 1);
	nbl_bord_cell.allocate(0, nti - 1, 0, ntj - 1);
	carte.allocate(0, nti - 1, 0, ntj - 1);// carte
	itopo.allocate(0, nti - 1, 0, ntj - 1);// carte
	carte.initialize();
	itopo.initialize(); itopo = 1; //default value

	//indices for regions
	reg_indices(param);
	
// lecture du fichier
	
	ifstream litcarte(param.str_file_mask.c_str(), ios::in);

	if (!litcarte) {cout<<endl<<"WARNING : Cannot read file: "<< param.str_file_mask.c_str() << endl;}
	for (int j=1; j<ntj-1; j++)
	{
		for (int i=1; i<nti-1; i++)
		{
			litcarte >> carte[i][j];
			if (carte[i][j]>nbl) carte[i][j] = nbl;

		}
	}
	litcarte.close();

	ifstream littopo(param.str_file_topo.c_str(), ios::in);

	if (!littopo) {cout<<endl<<"WARNING : Cannot read file: "<< param.str_file_topo.c_str() << endl;}
	for (int j=1; j<ntj-1; j++)
	{
		for (int i=1; i<nti-1; i++)
		{
			littopo >> itopo[i][j];
			//necessary to avoid Ha=0
			if (carte[i][j] && itopo(i,j)==0) itopo(i,j)=0.0001;
		}
	}
	littopo.close();

	if(param.nb_EEZ)
	{
		// mask contenant les codes des Zones exclusives
		maskEEZ.allocate(0, nti-1, 0, ntj-1);
		maskEEZ.initialize();

		ifstream rtxt(param.str_file_maskEEZ.c_str(), ios::in);
		if (!rtxt) {cout<<endl<<"WARNING : Cannot read file: "<< 
			param.str_file_maskEEZ.c_str() << endl;}

		for (int j=1;j<ntj-1;j++)
			for (int i=1;i<nti-1;i++)
			{
				rtxt >> maskEEZ[i][j];
			}
		rtxt.close(); 
	}

	if(param.mpa_simulation)
	{
		// mask contenant les codes des Zones exclusives
		maskMPA.allocate(0, nti-1, 0, ntj-1);
		maskMPA.initialize();

		ifstream rtxt(param.str_file_maskMPA.c_str(), ios::in);
		if (!rtxt) {cout<<endl<<"WARNING : Cannot read file: "<< param.str_file_maskMPA.c_str() << endl;}
		for (int j=1;j<ntj-1;j++)
			for (int i=1;i<nti-1;i++)
			{
				rtxt >> maskMPA[i][j];
			}
		rtxt.close(); 
	}


	//initialize
	global = 0;//Attn: for optimization no global until fixed!
	deltaX = param.deltaX;
	double reso_x = param.deltaX/60.0;
	cout << imax << " " << imin << " " << nti << endl;
	domain_type((nti-2)*reso_x);
//global = 0;//ATTN: put global to zero to be able to run with summy fisheries. Re-set when needed.
	definit_cell_bords(nti,ntj); // appel de la fonction 
	definit_lim_infsup(nti,ntj); // appel de la fonction
}

void PMap::reg_indices(CParam& param)
{
	//Regional indices
	const int nb_regs = param.nb_region;
	regimin.allocate(0,nb_regs-1);
	regjmin.allocate(0,nb_regs-1);
	regimax.allocate(0,nb_regs-1);
	regjmax.allocate(0,nb_regs-1);
	for (int a =0; a < param.nb_region; a++){

		double lonmin = param.area[a]->lgmin; double lonmax = param.area[a]->lgmax;
		double latmin = param.area[a]->ltmax; double latmax = param.area[a]->ltmin;
		regimin[a] = param.lontoi(lonmin); regimax[a] = param.lontoi(lonmax);
		regjmin[a] = param.lattoj(latmin); regjmax[a] = param.lattoj(latmax);

//		regimin[a] = (int) ( (param.area[a]->lgmin - param.longitudeMin)*60.0/param.deltaX +1);
//		regimax[a] = (int) ( (param.area[a]->lgmax - param.longitudeMin)*60.0/param.deltaX +1);
//		regjmax[a] = (int) ( (param.latitudeMax - param.area[a]->ltmin)*60.0/param.deltaY +1);
//		regjmin[a] = (int) ( (param.latitudeMax - param.area[a]->ltmax)*60.0/param.deltaY +1);

//cout << a+1 << "\t" << param.area[a]->lgmin << "\t" << param.area[a]->lgmax << "\t" << param.area[a]->ltmin << "\t" << param.area[a]->ltmax << endl;
//cout << a+1 << "\t" << regimin[a] << "\t" << regimax[a] << "\t" << regjmin[a] << "\t" << regjmax[a] << endl;
	}
}
/*
char PMap::get_bord_layer_x(const int i, const int j){

	bord_layer.b	= nbl_bord_cell[i][j];
	char pos	= bord_layer.cotex();
	return pos;	
}

char PMap::get_bord_layer_y(const int i, const int j){

	bord_layer.b	= nbl_bord_cell[i][j];
	char pos	= bord_layer.cotey();
	return pos;	
}
*/
void PMap::definit_cell_bords(const int nti, const int ntj) // définition des bords des cellules de la zone
{
  CBord bordures; //bordures, objet de la classe Cbords  

  int i, j;
  
  // Initialisation des cellules avec tous les cotes ouverts
  // SANS = sans frontiere, tous cotes ouverts
  // G_FERME = ferme a gauche
  // D_FERME = ferme a droite

  bordures.b = SANS;
  for (i = 1; i < nti-1; i++)
  {
    for (j = 1; j < ntj-1; j++)
    {
      bord_cell[i][j] = bordures.b;
      nbl_bord_cell[i][j] = bordures.b;
    }
  }  
  //***********************
  // Cadre de la grille 
  // Definition des cotés des cellules situées aux bords de la grille
  //***********************
  // bord horizontal haut de la grille, donc ferme a gauche
  i = 1;
  bordures.cote.y = SANS;
  bordures.cote.x = G_FERME;
//cout << "i=1 " << (int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;
  for (j = 1; j < ntj-1; j++)
  {
    bord_cell[i][j] = bordures.b;
    nbl_bord_cell[i][j] = bordures.b;
  }
  
  // bord horizontal bas de la grille donc ferme a droite
  i = nti-2;
  bordures.cote.y = SANS;
  bordures.cote.x = D_FERME;
//cout << "i=nbi-2 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " " << bordures.b << endl;
  for (j = 1; j < ntj-1; j++)
  {
    bord_cell[i][j] = bordures.b;
    nbl_bord_cell[i][j] = bordures.b;
  }

  // bord vertical gauche de la grille donc ferme a gauche
  j = 1;
  bordures.cote.x = SANS;
  bordures.cote.y = G_FERME;
//cout << "j=1 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;
  for (i = 1; i < nti-1; i++)
  {
    bord_cell[i][j] = bordures.b;
    nbl_bord_cell[i][j] = bordures.b;
  }

  // bord vertical droit de la grille donc ferme a droite
  j = ntj-2;
  bordures.cote.x = SANS;
  bordures.cote.y = D_FERME;
//cout << "j=nbj-2 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;
  for (i = 1; i < nti-1; i++)
  {
    bord_cell[i][j] = bordures.b;
    nbl_bord_cell[i][j] = bordures.b;
  }

  // finir les coins du cadre de la grille
  
  //coin en haut a gauche  
  bordures.cote.x = G_FERME;
  bordures.cote.y = G_FERME;
  bord_cell[1][1] = bordures.b;
  nbl_bord_cell[1][1] = bordures.b;
//cout << "i=1; j=1 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;
  
  //coin en haut a droite
  bordures.cote.x = G_FERME;
  bordures.cote.y = D_FERME;
  bord_cell[1][ntj-2]=bordures.b;
  nbl_bord_cell[1][ntj-2]=bordures.b;
//cout << "i=1; j=nbj-2 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;
  
  // coin en bas a gauche
  bordures.cote.x = D_FERME;
  bordures.cote.y = G_FERME;
  bord_cell[nti-2][1] = bordures.b;
  nbl_bord_cell[nti-2][1] = bordures.b;
//cout << "i=nbi-2; j=1 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl;

  // coin en bas a droite
  bordures.cote.x = D_FERME;
  bordures.cote.y = D_FERME;
  bord_cell[nti-2][ntj-2] = bordures.b; 
  nbl_bord_cell[nti-2][ntj-2] = bordures.b; 
//cout << "i=nbi-2; j=nbj-2 " <<(int)bordures.cote.x <<" " <<(int)bordures.cote.y << " "<< bordures.b << endl; 
//exit(1);
  //************************
  // Interieur de la grille
  // cas ou il y a des terres isolees 
  //************************
  for (i = 1; i < nti-1; i++)  {
    for (j = 1; j < ntj-1; j++){
      if (!carte[i][j]) {// si la valeur de la carte = 0 (donc terre)
		if (i>1){	// si cellule n'est pas sur le bord horizontal haut
		
			bordures.b = bord_cell[i-1][j];
			if (bordures.cotex() == SANS)
			{
			bordures.cote.x = D_FERME;
			bord_cell[i-1][j] = bordures.b;
			}
		}
		if (i<nti-2){ // si cellule n'est pas sur le bord horizontal bas
		
			bordures.b = bord_cell[i+1][j];
			if (bordures.cotex() == SANS)
			{
			bordures.cote.x = G_FERME;
			bord_cell[i+1][j] = bordures.b;
			}
		}
		if (j>1){ // si cellule n'est pas sur le bord vertical gauche
		
			bordures.b = bord_cell[i][j-1];
			if (bordures.cotey() == SANS)
			{
			bordures.cote.y = D_FERME;
			bord_cell[i][j-1] = bordures.b;
			}
		}
		if (j<ntj-2){ // si cellule n'est pas sur le bord vertical droit
		
			bordures.b = bord_cell[i][j+1];
			if (bordures.cotey() == SANS)
			{
		bordures.cote.y = G_FERME;
		bord_cell[i][j+1] = bordures.b;
		  	}
		}
		bordures.cote.x = TERRE;
		bordures.cote.y = TERRE;
		bord_cell[i][j] = bordures.b;
      }
	} // fin de boucle j
  } // fin de boucle i

  //set the boundaries for nbl-layer cells
  for (i = 1; i < nti-1; i++)  {
    for (j = 1; j < ntj-1; j++){
      if (carte[i][j]<nbl) { // if number of layers in the cell < total nbl
		if (i>1){	// si cellule n'est pas sur le bord horizontal haut
		
			  bordures.b = nbl_bord_cell[i-1][j];
			  if (bordures.cotex() == SANS){
			    bordures.cote.x = D_FERME;
			    nbl_bord_cell[i-1][j] = bordures.b;
			  }
		}
		if (i<nti-2){ // si cellule n'est pas sur le bord horizontal bas
		
			bordures.b = nbl_bord_cell[i+1][j];
			if (bordures.cotex() == SANS){
			  bordures.cote.x = G_FERME;
			  nbl_bord_cell[i+1][j] = bordures.b;
			}
		}
		if (j>1){ // si cellule n'est pas sur le bord vertical gauche
		
			bordures.b = nbl_bord_cell[i][j-1];
			if (bordures.cotey() == SANS){
			  bordures.cote.y = D_FERME;
			  nbl_bord_cell[i][j-1] = bordures.b;
			}
		}
		if (j<ntj-2){ // si cellule n'est pas sur le bord vertical droit
		
			bordures.b = nbl_bord_cell[i][j+1];
			if (bordures.cotey() == SANS){
	   		  bordures.cote.y = G_FERME;
			  nbl_bord_cell[i][j+1] = bordures.b;
		  	}
		}
		bordures.cote.x = TERRE;
		bordures.cote.y = TERRE;
		nbl_bord_cell[i][j] = bordures.b;
      }
	} // fin de boucle j
  } // fin de boucle i


  //******************
  // cas ou il ya des Terres dans les bords de la grille
  //******************
  i=1;
  		for (j=1 ; j < ntj-1 ; j++){
  			bordures.b=bord_cell[i][j];
  			if (carte[i+1][j]==0) {
	  			bordures.cote.x=TERRE;
	  			bord_cell[i][j]=bordures.b;
			}
		}
		
  i=nti-2;
  		for (j=1 ; j < ntj-1 ; j++){
  			bordures.b=bord_cell[i][j];
  			if (carte[i-1][j]==0) {
	  			bordures.cote.x=TERRE;
	  			bord_cell[i][j]=bordures.b;
	  		}
  		}

  j=1;
  		for (i=1 ; i < nti-1 ; i++){
  			bordures.b=bord_cell[i][j];
  			if (carte[i][j+1]==0) {
	  			bordures.cote.y=TERRE;
	  			bord_cell[i][j]=bordures.b;
	  		}
	  	}
	  	
  j=ntj-2;
  		for (i=1 ; i < nti-1 ; i++){
  			bordures.b=bord_cell[i][j];
  			if (carte[i][j-1]==0) {
	  			bordures.cote.y=TERRE;
	  			bord_cell[i][j]=bordures.b;
	  		}
	  	}

  //once more, similarly let's set boundaries for nb_layer cells
  i=1;
  		for (j=1 ; j < ntj-1 ; j++){
  			bordures.b=nbl_bord_cell[i][j];
  			if (carte[i+1][j]<nbl) {
	  			bordures.cote.x=TERRE;
	  			nbl_bord_cell[i][j]=bordures.b;
			}
		}
		
  i=nti-2;
  		for (j=1 ; j < ntj-1 ; j++){
  			bordures.b=nbl_bord_cell[i][j];
  			if (carte[i-1][j]<nbl) {
	  			bordures.cote.x=TERRE;
	  			nbl_bord_cell[i][j]=bordures.b;
	  		}
  		}
        //Inna Oct05: open left-right boundary
        if (global){
                i = 1;
                for (j=1 ; j < ntj-1 ; j++){
                        bordures.b=bord_cell[i][j];
                        if (carte[1][j] && carte[nti-2][j]){
                                bordures.cote.x=SANS;
                                bord_cell[i][j]=bordures.b;
                        }
                }

                i=nti-2;
                for (j=1 ; j < ntj-1 ; j++){
                        bordures.b=bord_cell[i][j];
                        if (carte[1][j] && carte[nti-2][j])
                        {
                                bordures.cote.x=SANS;
                                bord_cell[i][j]=bordures.b;
                        }
                }
        }
//END OF TEST		
		

  j=1;
  		for (i=1 ; i < nti-1 ; i++){
  			bordures.b=nbl_bord_cell[i][j];
  			if (carte[i][j+1]<nbl) {
	  			bordures.cote.y=TERRE;
	  			nbl_bord_cell[i][j]=bordures.b;
	  		}
	  	}

	  	
  j=ntj-2;
  		for (i=1 ; i < nti-1 ; i++){
  			bordures.b=nbl_bord_cell[i][j];
  			if (carte[i][j-1]<nbl) {
	  			bordures.cote.y=TERRE;
	  			nbl_bord_cell[i][j]=bordures.b;
	  		}
	  	}
}   


//***********************************************
// determine les valeurs d'indices a partir desquels 
// il n'y a plus de terre en partant des bords de la grille
// But: Optimization du calcul car inutile de calculer 
// le transport dans des cellules terre
//************************************************
void PMap::definit_lim_infsup(const int nti, const int ntj)

{
  const int rowmin = carte.rowmin();
  const int rowmax = carte.rowmax();
  const int colmin = carte.colmin();
  const int colmax = carte.colmax();

//cout << "rowmin " << rowmin << endl;
//cout << "rowmax " << rowmax << endl;
//cout << "colmin " << colmin << endl;
//cout << "colmax " << colmax << endl;

  IVECTOR iinftmp(colmin, colmax);  
  IVECTOR isuptmp(colmin, colmax);  
  iinftmp.initialize();
  isuptmp.initialize();

//cout << "iinftmp " << iinftmp << endl;
//cout << "isuptmp " << isuptmp << endl;

  for (int j = colmin; j <= colmax; j++)
  {
    int i = rowmin;
    while (i <= rowmax && carte[i][j] == 0)
      { i++; }
    iinftmp[j] = i;
  }

//cout << "iinftmp " << iinftmp << endl;

  for (int j = colmax; j >= colmin; j--)
  {
    int i = rowmax;
    while (i >= rowmin && carte[i][j] == 0)
      { i--; }
    isuptmp[j] = i;
  }

//cout << "isuptmp" << isuptmp << endl;

  jmin = colmin;
  while (iinftmp(jmin) > isuptmp(jmin)) jmin++;

  jmax = colmax;
  while (iinftmp(jmax) > isuptmp(jmax)) jmax--;

//cout << "jmin,jmax " << jmin << ' ' << jmax << endl;

 for (int j = jmin; j <= jmax; j++)
  {
	int i = iinftmp[j];
	if (!(carte[i - 1][j] == 0 && carte[i][j] != 0))
	{
		//cout << carte[i][j] << ' ' << carte[i + 1][j] << endl;
		//cout << i << ' ' << j << endl;
		//cout << __FILE__ << ':' << __LINE__ << endl;
		exit(1);
	}

	i = isuptmp[j];
	if (!(carte[i][j] != 0 && carte[i + 1][j] == 0))
	{
		//cout << carte[i][j] << ' ' << carte[i + 1][j] << endl;
		//cout << i << ' ' << j << endl;
		//cout << __FILE__ << ':' << __LINE__ << endl;
		exit(1);
	}

  }

  iinf.allocate(jmin, jmax);
  isup.allocate(jmin, jmax);
  iinf.initialize();
  isup.initialize();

  for (int j = jmin; j <= jmax; j++)
  {
	  iinf[j] = iinftmp[j];
	  isup[j] = isuptmp[j];
  }
//cout << "iinf " << iinf << endl;
//cout << "isup " << isup << endl;
//exit(1);
//=====================================

  IVECTOR jinftmp(rowmin, rowmax);  
  IVECTOR jsuptmp(rowmin, rowmax);  
  jinftmp.initialize();
  jsuptmp.initialize();
//cout << "jinftmp " << jinftmp << endl;
//cout << "jsuptmp " << jsuptmp << endl;

  for (int i = rowmin; i <= rowmax; i++)
  {
    int j = colmin;
    while (j <= colmax && carte[i][j] == 0)
      { j++; }
    jinftmp[i] = j;
  }
//cout << "jinftmp " << jinftmp << endl;
//exit(1);

  for (int i = rowmax; i >= rowmin; i--)
  {
    int j = colmax;
    while (j >= colmin && carte[i][j] == 0)
      { j--; }
    jsuptmp[i] = j;
  }
//cout << "jsuptmp " << jsuptmp << endl;
//exit(1);

  imin = rowmin;
  while (jinftmp(imin) > jsuptmp(imin)) imin++;
  imax = rowmax;
  while (jinftmp(imax) > jsuptmp(imax)) imax--;
//cout << imin << endl;
//cout << imax << endl;
//exit(1);

  jinf.allocate(imin, imax);
  jsup.allocate(imin, imax);
  for (int i = imin; i <= imax; i++)
  {
    jinf[i] = jinftmp[i];
    jsup[i] = jsuptmp[i];
  }
//cout << "jinf " << jinf << endl;
//cout << "jsup " << jsup << endl;

  for (int j=jmin; j <= jmax; j++)
  {
    int lb = iinf[j];
    int ub = isup[j];
    for (int i = lb; i <= ub; i++)
    {
      if (jinf[i] > j)
      {
        jinf[i] = j;
      }
      if (jsup[i] < j)
      {
	jsup[i] = j; 
      }
    }
  }

//cout << "jinf " << jinf << endl;
//cout << "jsup" << jsup << endl;
//exit(1);

  for (int i = imin; i <= imax; i++)
  {
    int lb = jinf[i];
    int ub = jsup[i];
    for (int j = lb; j <= ub; j++)
    {
      if (iinf[j] > i)
      {
        iinf[j] = i;
      }
      if (isup[j] < i)
      {
	isup[j] = i; 

      }
    }
  }

	// determining extended mask edges for matrices 
	// used in computing finite differences on bounds
  	imin1 = imin-1;
	int jmin1 = jmin-1;
	if (imin1 < rowmin) imin1 = rowmin;
	if (jmin1 < colmin) jmin1 = colmin;
	imax1 = imax+1;
	int jmax1 = jmax+1;
	if (imax1 > rowmax) imax1 = rowmax;
	if (jmax1 > colmax) jmax1 = colmax;

	jinf1.allocate(imin1,imax1); jinf1.initialize();
	jsup1.allocate(imin1,imax1); jsup1.initialize();
	
	for (int j = colmin; j <= colmax; j++){
		int i = rowmin;
		while (i < rowmax && carte[i+1][j] == 0) i++;
		iinftmp[j] = i;
	}

	for (int j = colmax; j >= colmin; j--){
		int i = rowmax;
		while (i > rowmin && carte[i-1][j] == 0) i--;
		isuptmp[j] = i;
	}

	for (int i = imin1; i <= imax1; i++){
		jinf1[i] = jinftmp[i];
		jsup1[i] = jsuptmp[i];
	}

	for (int j = jmin1; j <= jmax1; j++) {
		int lb = iinftmp[j];
		int ub = isuptmp[j];
		for (int i = lb; i <= ub; i++){
			if (jinf1[i] > j-1){
        			jinf1[i] = j-1;
      			}
      			if (jsup1[i] < j+1){
				jsup1[i] = j+1; 
      			}
    		}
  	}
	/////////////////////////////////////////////////

//cout << rowmin << " " << rowmax << " " << colmin << " " << colmax << endl;
//cout << imin << " " << imax << " " << imin1 << " " << imax1 << endl;
//cout << "jinf " << jinf << endl;
//cout << "jsup" << jsup << endl;
//cout << "jinf1 " << jinf1 << endl;
//cout << "jsup1" << jsup1 << endl;
//exit(1);

/*
  DMATRIX dmi(imin, imax, jinf, jsup);
  dmi.initialize();
  int count = 0;
  for (int i = imin; i <= imax; i++)
  {
    const int nonzeros = jsup[i] - jinf[i] + 1;
    count += jmax - nonzeros;
cout << i << ' ' << nonzeros << ' ' << count << endl;
    for (int j = jinf[i]; j <= jsup[i]; j++)
    {
      if (iinf[j] <= i && i <= isup[j])
      {
        dmi[i][j] += 1;
      }
    }
  }
//cout << dmi << endl;
cout << imin << endl;
cout << imax << endl;
cout << jmin << endl;
cout << jmax << endl;
cout << count << endl;
exit(1);
  DMATRIX dmj(jmin, jmax, iinf, isup);
  dmj.initialize();
  for (int j = jmin; j <= jmax; j++)
  {
    for (int i = iinf[j]; i <= isup[j]; i++)
    {
      if (jinf[i] <= j && j <= jsup[i])
      {
        dmj[j][i] += 1;
      }
    }
  }
cout << dmj << endl;
exit(1);
*/
}

void PMap::domain_type(const int nlon) {
	global = 0;
	if (nlon==360) 
		global = 1;

	cout << "Is global domain? ";
	if (global) cout << "YES" << endl;	
	else cout << "NO" << endl;	
}

