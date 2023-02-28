#ifndef __TRACE__
#define __TRACE__
#include <fstream>
#if (__GNUC__ >=3)
  using std::ofstream;
  using std::ostream;
#endif

#undef HERE
#define HERE cerr << "reached " << __LINE__ << " in " << __FILE__ << endl;

#undef TRACE
#define TRACE(object) cerr << "line " << __LINE__ << ", file " << __FILE__ << ", " << #object " = " << object << endl;

#undef mTRACE
#define mTRACE(object) cerr << "line " << __LINE__ << ", file " << __FILE__ << ", " << #object " = \n" << object << endl;

#undef TTRACE
#define TTRACE(o1,o2) cerr << "line " << __LINE__ << ", file " << __FILE__ << ", " << #o1 " = " << o1<< ", " << #o2 " = " << o2 << endl;

#undef TTTRACE
#define TTTRACE(o1,o2,o3) cerr << "line " << __LINE__ << ", file " << __FILE__ << ", " << #o1 " = " << o1<< ", " << #o2 " = " << o2 << ", " << #o3 " = " << o3 << endl;

#endif
