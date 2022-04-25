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
\brief Class has information about spatial domain: the land mask, the indexing and the boundaries.
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

// Classe permettant de gérer les cotés des cellules 
// dans les deux directions	x (i) et y (j)
// il ya 4 possibilites sur chaque axe: ouvert, ferme, ferme a droite, ferme a gauche
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
