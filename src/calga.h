#ifndef _calga_h_
#define _calga_h_
#include <cmath>
#include <cerrno>

double const defdokl = 1e-10;
//int lcalg5;

template < class T > double
calg5 (T & t, double (T::*f) (double), double x1, double x2, double dokl, int ile = 5)
{				//double static const sqrt2=sqrt(2);
  if (ile == 1)
    {
      static double const a = 0.53846931010568309105;	//sqrt(245-14*sqrt(70))/21
      static double const b = 0.90617984593866399282;	//sqrt(245+14*sqrt(70))/21
      static double const f0 = 128. / 225.;
      static double const fb = (1. / 3. - a * a * (1. - f0 / 2.)) / (b * b - a * a);	// 0.23692689
      static double const fa = 1. - f0 / 2. - fb;	// 0.47862867
      double h = (x2 - x1) / 2;
//      if (dokl == defdokl)
//	lcalg5 = 0;
      double xs = x1 + h;
      double y1m = (t.*f) (xs - h * a);
      double y1p = (t.*f) (xs + h * a);
      double y2m = (t.*f) (xs - h * b);
      double y2p = (t.*f) (xs + h * b);
      double ys = (t.*f) (xs);
      double r1 = (y1m + y1p) / 2 - ys;
      double r2 = (y2m + y2p) / 2 - ys;
//      lcalg5 += 5;
      if (fabs ((r2 * a * a - r1 * b * b) * 10 * h) < dokl)
	return h * (f0 * ys + fa * (y1m + y1p) + fb * (y2m + y2p));
      else
	{
	  dokl /= 2;
	  h /= 2;
	  return calg5 (t, f, x1, x1 + h, dokl) 
	       + calg5 (t, f, x1 + h, xs, dokl) 
	       + calg5 (t, f, xs, xs + h, dokl)
  	       + calg5 (t, f, xs + h, x2, dokl);
	}
    }
  else
    {
      double d = (x2 - x1) / ile;
      dokl /= sqrt (double(ile));
      double suma = 0;
      while (ile-- > 0)
	{
	  suma += calg5 (t, f, x1 + ile * d, x1 + ile * d + d, dokl, 1);
	}
      return suma;
    }
}

template < class T > double
calg5dokl (T & t, double (T::*f) (double, double), double x1, double x2,
	   double dokl = defdokl, int ile = 1)
{				//double static const sqrt2=sqrt(2);
  if (ile == 1)
    {
      static double const a = 0.53846931010568309105;	//sqrt(245-14*sqrt(70))/21
      static double const b = 0.90617984593866399282;	//sqrt(245+14*sqrt(70))/21
      static double const f0 = 128.0 / 225;
      static double const fb = (1.0 / 3 - a * a * (1 - f0 / 2)) / (b * b - a * a);	// 0.23692689
      static double const fa = 1 - f0 / 2 - fb;	// 0.47862867
      double h = (x2 - x1) / 2;
//      if (dokl == defdokl)
//	lcalg5 = 0;
      double xs = x1 + h;
      double y1m = (t.*f) (xs - h * a, dokl / h / 5);
      double y1p = (t.*f) (xs + h * a, dokl / h / 5);
      double y2m = (t.*f) (xs - h * b, dokl / h / 5);
      double y2p = (t.*f) (xs + h * b, dokl / h / 5);
      double ys = (t.*f) (xs, dokl / h / 5);
      double r1 = (y1m + y1p) / 2 - ys;
      double r2 = (y2m + y2p) / 2 - ys;
//      lcalg5 += 5;
      if (fabs ((r2 * a * a - r1 * b * b) * 10 * h) < dokl)
	return h * (f0 * ys + fa * (y1m + y1p) + fb * (y2m + y2p));
      else
	{
	  dokl /= 2;
	  h /= 2;
	  return calg5dokl (t, f, x1, x1 + h, dokl) 
	       + calg5dokl (t, f, x1 + h, xs, dokl) 
	       + calg5dokl (t, f, xs, xs + h, dokl) 
	       + calg5dokl (t, f, xs + h, x2, dokl);
	}
    }
  else
    {
      double d = (x2 - x1) / ile;
      dokl /= sqrt (ile);
      double suma = 0;
      while (ile-- > 0)
	{
	  suma += calg5dokl (t, f, x1 + ile * d, x1 + ile * d + d, dokl, 1);
	}
      return suma;

    }
}


#endif
