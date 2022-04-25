#include "SeapodymDocConsole.h"

void SeapodymDocConsole::UpdateDisplay()
{
       if (t_count>1){
		cout<<t_count<<"\t| "<<date_str<<" |\t"<<sumP<<"\t|\t"<<sumFprime[0]<<"\t|\t"<<sumF[0] << endl;
       } else{
		cout << endl << "Building forage populations:" << endl;
                cout << "step\t|\tdate\t|\tPrim. Prod\t|\tForage Prod\t|\tForage Biomass"<< endl;
                cout << "----------------------------------------------------------------------------------------"<< endl;
       }
}
