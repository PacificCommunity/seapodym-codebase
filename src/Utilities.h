#ifndef __Utilities_h__
#define __Utilities_h__

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream> 
#include <cmath>
// #include <climits>
#include <cstdlib>
#include "ctrace.h"

using  namespace std;

class Utilities
{
public:
	static inline string MakeDate(int yr, int mo, int jr)
	{
		ostringstream osyr;
		ostringstream osjr;
		osyr << yr;
		osjr << jr;
		string date= MonthName(mo)+" "+osjr.str()+", "+osyr.str()+"  ";
		return date;
	}

	static inline string MakeDate(int yr, int mo)
	{
		ostringstream osyr;
		osyr << yr;
		string date= MonthName(mo)+", "+osyr.str();
		return date;
	}

	static inline string MonthName(int mo)
	{
		switch (mo)
		{
		case 1:
			return "Jan";
			break;
		case 2:
			return "Feb";
			break;
		case 3:
			return "Mar";
			break;
		case 4:
			return "Apr";
			break;
		case 5:
			return "May";
			break;
		case 6:
			return "Jun";
			break;
		case 7:
			return "Jul";
			break;
		case 8:
			return "Aug";
			break;
		case 9:
			return "Sep";
			break;
		case 10:
			return "Oct";
			break;
		case 11:
			return "Nov";
			break;
		case 12:
			return "Dec";
			break;
		}
return "";
	}

// Converts an integer to a string.
	static inline string itoa(int i)
	{
		// Create an output stream.
		ostringstream oss;
		// Write into the stream.
		oss << i;
		// Gets the result string;
		return oss.str();
	}

	static inline int MyMax(int a, int b)
	{
		return (a > b)? a : b;
	}

	static inline double MyMax(double a, double b)
	{
		return (a > b)? a : b;
	}

	static inline short MyMax(short a, short b)
	{
		return (a > b)? a : b;
	}

	static inline char MyMax(char a, char b)
	{
		return (a > b)? a : b;
	}

    static inline int MyMin(int a, int b)
	{
		return (a < b)? a : b;
	}

    static inline double MyMin(double a, double b)
	{
		return (a < b)? a : b;
	}

   static inline short MyMin(short a, short b)
   {
		return (a < b)? a : b;
	}

    static inline char MyMin(char a, char b)
	{
		return (a < b)? a : b;
	}

   static inline int* create1d(int *mat, const int n1,  const int val=0)
	{
		mat=new int [n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=val;
		}
		return mat;
	}
	//create array of size (n1,n2,...) and initialise it with value val
   static inline double* create1d(double *mat,int n1, double val=0)
	{
		mat=new double[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=val;
		}
		return mat;
	}

	//create array of size (n1,n2,...) and initialise it with value val
   static inline string* create1d(string *mat,int n1, string val="")
	{
		mat=new string[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=val;
		}
		return mat;
	}

    static inline double** create2d(double **mat,int n1, int n2, double val=0)
	{
		mat=new double*[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double[n2];
			for (int j=0 ; j<n2 ; j++)
			{
				mat[i][j]=val;
			}
		}
		return mat;
	}


    static inline double** create2d(double **mat,int n1, const IVECTOR& n2, double val=0)
	{
		mat=new double*[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=val;
			}
		}
		return mat;
	}

    static inline string** create2d(string **mat,int n1, const IVECTOR& n2, string val="")
	{
		mat=new string*[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new string[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=val;
			}
		}
		return mat;
	}

    static inline int** create2d(int **mat,int n1, int n2, int val=0)
	{
		mat=new int*[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new int[n2];
			for (int j=0 ; j<n2 ; j++)
			{
				mat[i][j]=val;
			}
		}
		return mat;
	}

    static inline int** create2d(int **mat,int n1, const IVECTOR& n2, int val=0)
	{
		mat=new int*[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new int[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=val;
			}
		}
		return mat;
	}

    static inline double*** create3d(double ***mat,int n1, int n2, int n3, double val=0)
	{
		mat=new double**[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double*[n2];
			for (int j=0 ; j<n2 ; j++)
			{
				mat[i][j]=new double[n3];
				for (int k=0 ; k<n3 ; k++)
				{
					mat[i][j][k]=val;
				}
			}
		}
		return mat;
	}

   static inline double*** create3d(double ***mat,int n1, const IVECTOR& n2, const IVECTOR& n3, double val=0)
	{
		mat=new double**[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double*[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=new double[n3[i]];
				for (int k=0 ; k<n3[i] ; k++)
				{
					mat[i][j][k]=val;
				}
			}
		}
		return mat;
	}

    static inline int*** create3d(int ***mat,int n1, int n2, int n3, int val=0)
	{
		mat=new int**[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new int*[n2];
			for (int j=0 ; j<n2 ; j++)
			{
				mat[i][j]=new int[n3];
				for (int k=0 ; k<n3 ; k++)
				{
					mat[i][j][k]=val;
				}
			}
		}
		return mat;
	}

    static inline double**** create4d(double ****mat,int n1, int n2, int n3, int n4, double val=0)
	{
		mat=new double***[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double**[n2];
			for (int j=0 ; j<n2 ; j++)
			{
				mat[i][j]=new double*[n3];
				for (int k=0 ; k<n3 ; k++)
				{
					mat[i][j][k]=new double[n4];
					for (int l=0 ; l<n4 ; l++)
					mat[i][j][k][l]=val;
				}
			}
		}
		return mat;
	}

   static inline double**** create4d(double ****mat,int n1, const IVECTOR& n2, int n3, int n4, double val=0)
	{
		mat=new double***[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double**[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=new double*[n3];
				for (int k=0 ; k<n3 ; k++)
				{
					mat[i][j][k]=new double[n4];
					for (int l=0 ; l<n4 ; l++)
					mat[i][j][k][l]=val;
				}
			}
		}
		return mat;
	}

   static inline double**** create4d(double ****mat,int n1, const IVECTOR& n2, const IVECTOR& n3, const IVECTOR& n4, double val=0)
	{
		mat=new double***[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double**[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=new double*[n3[i]];
				for (int k=0 ; k<n3[i] ; k++)
				{
					mat[i][j][k]=new double[n4[i]];
					for (int l=0 ; l<n4[i] ; l++)
					mat[i][j][k][l]=val;
				}
			}
		}
		return mat;
	}

   static inline double***** create5d(double *****mat,int n1, const IVECTOR& n2, int n3, const IVECTOR& n4, const IVECTOR& n5, double val=0)
	{
		mat=new double****[n1];
		for (int i=0 ; i<n1 ; i++)
		{
			mat[i]=new double***[n2[i]];
			for (int j=0 ; j<n2[i] ; j++)
			{
				mat[i][j]=new double**[n3];
				for (int k=0 ; k<n3 ; k++)
				{
					mat[i][j][k]=new double*[n4[i]+1];
					for (int l=0 ; l<(n4[i]+1) ; l++)
					{
						mat[i][j][k][l]=new double[n5[i]+1];
						for (int m=0 ; m<(n5[i]+1) ; m++)
						{
							mat[i][j][k][l][m]=val;
						}
					}
				}
			}
		}
		return mat;
	}

    static inline void delete1d(string *mat)
	{
		delete [] mat;
		mat = 0;
	}
    static inline void delete1d(const IVECTOR& mat)
	{
	}
    static inline void delete1d(double *mat)
	{
		delete [] mat;
		mat = 0;
	}

	//delete array of size (n1,n2,...)
    static inline void delete2d(double **mat,int n1)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete2d(int **mat,int n1)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete2d(string **mat,int n1)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			delete [] mat[i];
		}
		delete [] mat;
	}


    static inline void delete3d(double ***mat,int n1,int n2)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2 ; j++)
			{
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete3d(double ***mat,int n1,const IVECTOR& n2)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2[i] ; j++)
			{
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete3d(int ***mat,int n1,int n2)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2 ; j++)
			{
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete4d(double ****mat,int n1,int n2,int n3)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2 ; j++)
			{
				for (int k=0 ; k<n3 ; k++)
				{
					delete [] mat[i][j][k];
				}
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete4d(double ****mat,int n1,const IVECTOR& n2,int n3)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2[i] ; j++)
			{
				for (int k=0 ; k<n3 ; k++)
				{
					delete [] mat[i][j][k];
				}
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete4d(double ****mat,int n1,const IVECTOR& n2,const IVECTOR& n3)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2[i] ; j++)
			{
				for (int k=0 ; k<n3[i] ; k++)
				{
					delete [] mat[i][j][k];
				}
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

    static inline void delete5d(double *****mat,int n1,const IVECTOR& n2,int n3,const IVECTOR& n4)
	{
		for (int i=0 ; i<n1 ; i++)
		{
			for (int j=0 ; j<n2[i] ; j++)
			{
				for (int k=0 ; k<n3 ; k++)
				{
					for (int l=0 ; l<(n4[i]+1) ; l++)
					{
						delete [] mat[i][j][k][l];
					}
					delete [] mat[i][j][k];
				}
				delete [] mat[i][j];
			}
			delete [] mat[i];
		}
		delete [] mat;
	}

};
#endif
