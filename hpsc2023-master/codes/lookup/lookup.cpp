#include <iostream>
using std :: cout;
using std :: endl;

double lookupVal(int n, double *x, double *y, double xval)
{
  for ( int i = 0 ; i < n ; ++i)
      if ( xval >= x[i] && xval <= x[i+1] )
	  return y[i] + (xval - x[i]) * (y[i+1]-y[i]) / (x[i+1]-x[i]);
  return 0.;
}

int main(int argc, char *argv[])
{

  int numPE, myPE;

  int n = 100;
  double x[n], y[n];

  for ( int i = 0 ; i < n ; ++i)
    {
      x[i] = i;
      y[i] = i*i;
    }

  double xval = 2.5;  double yval = lookupVal(n,x,y,xval);

  cout << "For xval = " << xval << " yval = " << yval << endl;

  return 0;
  
}
