#ifndef __MAP_H__
#define __MAP_H__

#include "Param.h"

enum type_de_bords
{
	TERRE=-1,
	SANS=0,
	G_FERME,
	D_FERME
};

/*!
\brief Class managing spatial domain and grid: the land mask, the indexing and the boundaries.
\details This class reads land mask, EEZ mask (if exist) and topographic indices. 
The boundary conditions are defined here as well using the land mask information.
Also, the ragged array indices are computed and stored in this class.
*/
class PMap  
{
public:
	PMap() {/*DoesNothing*/};
	virtual ~PMap() {/*DoesNothing*/};

public:
	void lit_map(CParam &param);
	void delete_map(const CParam &param) {/*DoesNothing*/}
	void reg_indices(CParam &param);

	//char get_bord_layer_x(const int i, const int j);
	//char get_bord_layer_y(const int i, const int j);

public:

	IMATRIX bord_cell;	// bord des cellules (voir Classe CBord plus loin)
	IMATRIX nbl_bord_cell;	// border for the cells with all layers being exist
	IMATRIX carte;		// carte (mettre des zéros quand terre ou ile)
	DMATRIX itopo;		// carte (mettre des zéros quand terre ou ile)
	IMATRIX maskEEZ;	// carte des EEZ avec un code par EEZ
	IMATRIX maskMPA;	// carte des MPA avec un code par MPA


	int imin;
	int imax;
	int jmin;
	int jmax;
	int imin1;
	int imax1;
	int global;

	IVECTOR iinf;
	IVECTOR isup;
	IVECTOR jinf;
	IVECTOR jsup;
	IVECTOR jinf1;
	IVECTOR jsup1;

	ivector regimin;
	ivector regimax;
	ivector regjmin;
	ivector regjmax;

private:
	void definit_cell_bords(const int nti, const int ntj);	// définition des bords des cellules
	void definit_lim_infsup(const int nti, const int ntj);	// définition des limites inférieures 
								// et supérieures de la matrice zone pour optimiser le calcul
	void domain_type(const int nlon); 
};

/*!
\brief Class handling the type of the borders of a grid cell.
\details For a given pair of indices (i,j) structure cote stores the type of cell's borders, two in x and two in y direction - left-closed (G_FERME), right-closed (D_FERME) or open (SANS) for ocean cells, and land (TERRE) if the land is next to the land cell. 
*/
class CBord
{
public:
	union
	{
		unsigned short int b;
		struct
		{
			char x;
			char y;
		}cote;
	};

	int cotex() {return (int)cote.x;}
	int cotey() {return (int)cote.y;}
};
#endif //  __MAP_H__
