#ifndef __Cmdopt_h__
#define __Cmdopt_h__

#include <iostream>

using namespace std;

/*!
\Class written to handle command line arguments/
*/
class Cmdopt
{
    public:
        Cmdopt(int argc, char** argv);
        virtual ~Cmdopt()  {/*DoNothing*/};

    public:
        int cmp_regime = -1;
        int sftype = 0;
        int nb_aat = 1;
        bool reset_buffers = false;
        bool file_exists = false;

        void help(char* argv0);
        void help_sim(char* argv0);
        void CheckInfo(char* Option, char* argv0);
        void OptionToCode(char* Option);
        void read_sa_options(int argc, char** argv);
        char* getCmdOption(char ** begin, char ** end, const std::string& option);
        bool cmdOptionExists(char** begin, char** end, const std::string& option);
};

#endif