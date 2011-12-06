//by JN
#include <cmath>
#include <cstdlib>

int
charge1 (int code1)
{
  double u, c, t;
  double d, s, b;
  double q1c = 0, q2c = 0, q3c = 0;
  double sum;
  int code = abs (code1);

  u = c = t = 2;		//2,4,6
  d = s = b = -1;		//1,3,5

  int q1 = (code / 1000) % 10;
  int q2 = (code / 100) % 10;
  int q3 = (code / 10) % 10;

  if (q1 == 0)
    {
      q1c = 0;
    }
  if (q1 == 2 || q1 == 4 || q1 == 6)
    {
      q1c = u;
    }
  if (q1 == 1 || q1 == 3 || q1 == 5)
    {
      q1c = d;
    }

  if (q2 == 0)
    {
      q2c = 0;
    }
  if (q2 == 2 || q2 == 4 || q2 == 6)
    {
      q2c = u;
    }
  if (q2 == 1 || q2 == 3 || q2 == 5)
    {
      q2c = d;
    }

  if (q3 == 0)
    {
      q3c = 0;
    }
  if (q3 == 2 || q3 == 4 || q3 == 6)
    {
      q3c = u;
    }
  if (q3 == 1 || q3 == 3 || q3 == 5)
    {
      q3c = d;
    }

  if (q1 != 0)
    {
      sum = (q1c + q2c + q3c) / 3.;
    }
  else
    {
      sum = (-1 * q2c + q3c) / 3.;
      if (sum < 0)
	{
	  sum = -1 * sum;
	}
    }

//cout<<q1<<" "<<q2<<" "<<q3<<endl;
//cout<<q1c<<" "<<q2c<<" "<<q3c<<endl;

  int charge = int (sum);

  if (code1 > 0)
    {
      return charge;
    }
  else
    {
      return -1 * charge;
    }
}
