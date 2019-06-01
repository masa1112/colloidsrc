#ifndef RAN
#define RAN
#define IA 16807
#define IM 2147483647
//#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define EPS (1.2E-07)
#define MAX(a,b) (a>b)?a:b
#define MIN(a,b) (a<b)?a:b

#include<time.h>
/*
inline double ran()
{
  static int initial=0;
  static long idum;
  time_t m_time;
  struct tm *jikoku;
  long k;
  double ans;

  if(initial==0){
    initial=1;
    time(&m_time);
    jikoku=localtime(&m_time);
    idum=jikoku[0].tm_sec;
    //printf("ran begins.\n");
  }    
  
  idum ^= MASK;
  k=(idum)/IQ;
  idum=IA*(idum-k*IQ)-IR*k;
  if (idum < 0) idum += IM;
  ans=AM*(idum);
  idum ^= MASK;
  return ans;
}

#endif*/

/* (C) Copr. 1986-92 Numerical Recipes Software t+!-). */

/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */


inline double ran()
//int *idum;
{
  static int initial=0;
  static long idum;
  time_t m_time;
  struct tm *jikoku;
  //long k;
  //double ans;
  
  if(initial==0){
    initial=1;
    idum=(unsigned)time(NULL);
  }
  
  
  int j,k;
  static int iv[NTAB],iy=0;
  void nrerror();
  static double NDIV = 1.0/(1.0+(IM-1.0)/NTAB);
  static double RNMX = (1.0-EPS);
  static double AM = (1.0/IM);

	if ((idum <= 0) || (iy == 0)) {
		idum = MAX(-idum,idum);
                for(j=NTAB+7;j>=0;j--) {
			k = idum/IQ;
			idum = IA*(idum-k*IQ)-IR*k;
			if(idum < 0) idum += IM;
			if(j < NTAB) iv[j] = idum;
		}
		iy = iv[0];
	}
	k = idum/IQ;
	idum = IA*(idum-k*IQ)-IR*k;
	if(idum<0) idum += IM;
	j = iy*NDIV;
	iy = iv[j];
	iv[j] = idum;
	return MIN(AM*iy,RNMX);
}
#undef IA 
#undef IM 
#undef IQ
#undef IR
#undef NTAB
#undef EPS 
#undef MAX
#undef MIN
#endif
