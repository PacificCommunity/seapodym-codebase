#include "Dimensions.h"
#include "XMLDocument2.h"


Dimensions::Dimensions()
{
	try {
		doc.read("parfile_new.xml");
	} catch (runtime_error e) {
		cerr << e.what() << '\n';
		//return false;
	}

		nbi = doc.getInteger("/deltaX", "value");;//10; 
		nbj = 0; nbz = 0; 
		//temporal 
		nbt = 0;
		//populations
		nb_forage = 0; nb_species = 0; nb_cohorts = 0;
		//cohorts
		nbc_lv = 0; nbc_jv = 0; nbc_yn = 0; nbc_ad = 0; i0_ad = 0;

}

Dimensions& Dimensions::getInstance()
{
  static Dimensions instance;
  return instance;
}





