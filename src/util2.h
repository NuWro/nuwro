#ifndef _util_h_
#define _util_h_ 
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <class F, class Y>
class Fix2nd
{ F _f;
  Y _y;
  public:
  Fix2nd(F  f, Y y) :_f(f),_y(y){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y);
  }
} ;

template <class F, class Y>
inline Fix2nd<F,Y> fix2nd(F  f, Y y) 
{
  return Fix2nd<F,Y>(f,y);
}

template <class F, class Y>
class Fix2
{ F _f;
  Y _y;
  public:
  Fix2(F  f, Y y) :_f(f),_y(y){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y);
  }

} ;

template <class F, class Y>
inline Fix2<F,Y> fix2(F  f, Y y) 
{
  return Fix2<F,Y>(f,y);
}


template <class F, class X>
class Fix1st
{ F _f;
  X _x;
  public:
  Fix1st(F  f, X x) :_f(f),_x(x){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y);
  }
} ;

template <class F, class X>
inline Fix1st<F,X> fix1st(F  f, X x) 
{
  return Fix1st<F,X>(f,x);
}


template <class F, class X>
class Fix1
{ F _f;
  X _x;
  public:
  Fix1(F  f, X x) :_f(f),_x(x){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y);
  }
};

template <class F, class X>
inline Fix1<F,X> fix1(F  f, X x) 
{
  return Fix1<F,X>(f,x);
}


template <class F, class X, class Y>
class Fix12
{ F _f;
  X _x;
  Y _y;
  public:
  Fix12(F  f, X x,Y y) :_f(f),_x(x),_y(y){}
  template <class Z>
  double operator()(Z z)
  { return _f(_x,_y,z);
  }
};

template <class F, class X, class Y>
inline Fix12<F,X,Y> fix12(F  f, X x, Y y) 
{
  return Fix12<F,X,Y>(f,x,y);
}

template <class F, class X, class Z>
class Fix13
{ F _f;
  X _x;
  Z _z;
  public:
  Fix13(F  f, X x,Z z) :_f(f),_x(x),_z(z){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y,_z);
  }

} ;

template <class F, class X, class Z>
inline Fix13<F,X,Z> fix13(F  f, X x, Z z) 
{
  return Fix13<F,X,Z>(f,x,z);
}

template <class F, class Y, class Z>
class Fix23
{ F _f;
  Y _y;
  Z _z;
  public:
  Fix23(F  f, Y y, Z z) :_f(f),_y(y),_z(z){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y,_z);
  }

} ;

template <class F, class Y, class Z>
inline Fix23<F,Y,Z> fix23(F  f, Y y, Z z) 
{
  return Fix23<F,Y,Z>(f,y,z);
}

template <class F,class X, class Y, class Z>
class Fix123
{ F _f;
  X _x;
  Y _y;
  Z _z;
  public:
  Fix123(F  f, X x, Y y, Z z) :_f(f),_x(x),_y(y),_z(z){}
  template <class T>
  double operator()(T t)
  { return _f(_x,_y,_z,t);
  }

} ;

template <class F, class X, class Y, class Z>
inline Fix123<F,X,Y,Z> fix123(F  f, X x, Y y, Z z) 
{
  return Fix123<F,X,Y,Z>(f,x,y,z);
}

template <class F, class X, class Y, class T>
class Fix124
{ F _f;
  X _x;
  Y _y;
  T _t;
  public:
  Fix124(F  f, X x, Y y, T t) :_f(f),_x(x),_y(y),_t(t){}
  template <class Z>
  double operator()(Z z)
  { return _f(_x,_y,z,_t);
  }

} ;

template <class F, class X, class Y, class T>
inline Fix124<F,X,Y,T> fix124(F  f, X x, Y y, T t) 
{
  return Fix124<F,X,Y,T>(f,x,y,t);
}

template <class F, class X, class Z, class T>
class Fix134
{ F _f;
  X _x;
  Z _z;
  T _t;
  public:
  Fix134(F f, X x, Z z, T t) :_f(f),_x(x),_z(z),_t(t){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y,_z,_t);
  }

} ;

template <class F, class X, class Z, class T>
inline Fix134<F,X,Z,T> fix134(F  f, X x, Z z, T t) 
{
  return Fix134<F,X,Z,T>(f,x,z,t);
}

template <class F, class Y, class Z, class T>
class Fix234
{ F _f;
  Y _y;
  Z _z;
  T _t;
  public:
  Fix234(F  f, Y y, Z z, T t) :_f(f),_y(y),_z(z),_t(t){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y,_z,_t);
  }

} ;

template <class F, class Y, class Z, class T>
inline Fix234<F,Y,Z,T> fix234(F  f, Y y, Z z, T t) 
{
  return Fix234<F,Y,Z,T>(f,y,z,t);
}


template <class F, class X, class Y, class Z, class T>
class Fix1234
{ F _f;
  X _x;
  Y _y;
  Z _z;
  T _t;
  public:
  Fix1234(F  f, X x, Y y, Z z, T t) :_f(f),_x(x),_y(y),_z(z),_t(t){}
  template <class U>
  double operator()(U u)
  { return _f(_x,_y,_z,_t,u);
  }
} ;

template <class F, class X, class Y, class Z, class T>
inline Fix1234<F,X,Y,Z,T> fix1234(F  f, X x, Y y, Z z, T t) 
{
  return Fix1234<F,X,Y,Z,T>(f,x,y,z,t);
}

template <class F, class X, class Y, class Z, class U>
class Fix1235
{ F _f;
  X _x;
  Y _y;
  Z _z;
  U _u;
  public:
  Fix1235(F  f, X x, Y y, Z z, U u) :_f(f),_x(x),_y(y),_z(z),_u(u){}
  template <class T>
  double operator()(T t)
  { return _f(_x,_y,_z,t,_u);
  }
} ;

template <class F, class X, class Y, class Z, class U>
inline Fix1235<F,X,Y,Z,U> fix1235(F  f, X x, Y y, Z z, U u) 
{
  return Fix1235<F,X,Y,Z,U>(f,x,y,z,u);
}

template <class F, class X, class Y, class T, class U>
class Fix1245
{ F _f;
  X _x;
  Y _y;
  T _t;
  U _u;
  public:
  Fix1245(F  f, X x, Y y, T t, U u) :_f(f),_x(x),_y(y),_t(t),_u(u){}
  template <class Z>
  double operator()(Z z)
  { return _f(_x,_y,z,_t,_u);
  }
} ;

template <class F, class X, class Y, class T, class U>
inline Fix1245<F,X,Y,T,U> fix1245(F  f, X x, Y y, T t, U u) 
{
  return Fix1245<F,X,Y,T,U>(f,x,y,t,u);
}




template <class F,  class X, class Z, class T, class U>
class Fix1345
{ F _f;
  X _x;
  Z _z;
  T _t;
  U _u;
  public:
  Fix1345(F  f, X x, Z z, T t, U u) :_f(f),_x(x),_z(z),_t(t),_u(u){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y,_z,_t,_u);
  }
} ;

template <class F , class X, class Z, class T, class U>
inline Fix1345<F,X,Z,T,U> fix1345(F  f, X x, Z z, T t, U u) 
{
  return Fix1345<F,X,Z,T,U>(f,x,z,t,u);
}

template <class F, class Y, class Z, class T, class U>
class Fix2345
{ F _f;
  Y _y;
  Z _z;
  T _t;
  U _u;
  public:
  Fix2345(F  f, Y y, Z z, T t, U u) :_f(f),_y(y),_z(z),_t(t),_u(u){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y,_z,_t,_u);
  }
} ;

template <class F, class Y, class Z, class T, class U>
inline Fix2345<F,Y,Z,T,U> fix2345(F  f, Y y, Z z, T t, U u) 
{
  return Fix2345<F,Y,Z,T,U>(f,y,z,t,u);
}



template <class _T,class F>
class Makefun
{ _T & _t;
  F _f;
 public: 
 Makefun(_T& t, F f):_t(t),_f(f){}
  template <class X>
 double operator()(X x)
 {return (_t.*_f)(x);
 }
  template <class X,class Y>
 double operator()(X x, Y y)
 {return (_t.*_f)(x,y);
 }
  template <class X,class Y,class Z>
double operator()(X x, Y y,Z z)
 {return (_t.*_f)(x,y,z);
 }
  template <class X,class Y,class Z, class T >
double operator()(X x, Y y,Z z,T t)
 {return (_t.*_f)(x,y,z,t);
 }
  template <class X,class Y,class Z, class T, class U >
double operator()(X x, Y y,Z z,T t,U u)
 {return (_t.*_f)(x,y,z,t,u);
 }
  template <class X,class Y,class Z, class T, class U,class V >
double operator()(X x, Y y,Z z,T t,U u,V v)
 {return (_t.*_f)(x,y,z,t,u,v);
 }
};

template <class T,class F>
Makefun<T,F>  makefun(T& t, F f)
 {return Makefun<T,F>(t,f);
 }

template <class F, class Y, class Z, class T, class U, class W>
class Fix23456
{ F _f;
  Y _y;
  Z _z;
  T _t;
  U _u;
  W _w;
  public:
  Fix23456(F  f, Y y, Z z, T t, U u, W w) :_f(f),_y(y),_z(z),_t(t),_u(u), _w(w){}
  template <class X>
  double operator()(X x)
  { return _f(x,_y,_z,_t,_u, _w);
  }
} ;


template <class F, class Y, class Z, class T, class U, class W>
inline Fix23456<F,Y,Z,T,U,W> fix23456(F  f, Y y, Z z, T t, U u, W w) 
{
  return Fix23456<F,Y,Z,T,U,W>(f,y,z,t,u,w);
}


template <class F, class X, class Z, class T, class U, class W>
class Fix13456
{ F _f;
  X _x;
  Z _z;
  T _t;
  U _u;
  W _w;
  public:
  Fix13456(F  f, X x, Z z, T t, U u, W w) :_f(f),_x(x),_z(z),_t(t),_u(u), _w(w){}
  template <class Y>
  double operator()(Y y)
  { return _f(_x,y,_z,_t,_u, _w);
  }
} ;



template <class F, class X, class Z, class T, class U, class W>
inline Fix13456<F,X,Z,T,U,W> fix13456(F  f, X x, Z z, T t, U u, W w) 
{
  return Fix13456<F,X,Z,T,U,W>(f,x,z,t,u,w);
}


template <class F, class X, class Y, class T, class U, class W>
class Fix12456
{ F _f;
  X _x;
  Y _y;
  T _t;
  U _u;
  W _w;
  public:
  Fix12456(F  f, X x, Y y, T t, U u, W w) :_f(f),_x(x),_y(y),_t(t),_u(u), _w(w){}
  template <class Z >
  double operator()(Z z)
  { return _f(_x,_y,z,_t,_u, _w);
  }
} ;



template <class F, class X, class Y, class T, class U, class W>
inline Fix12456<F,X,Y,T,U,W> fix12456(F  f, X x, Y y, T t, U u, W w) 
{
  return Fix12456<F,X,Y,T,U,W>(f,x,y,t,u,w);
}

template <class F, class X, class Y, class Z, class U, class W>
class Fix12356
{ F _f;
  X _x;
  Y _y;
  Z _z;
  U _u;
  W _w;
  public:
  Fix12356(F  f, X x, Y y, Z z, U u, W w) :_f(f),_x(x),_y(y),_z(z),_u(u), _w(w){}
  template <class T>
  double operator()(T t)
  { return _f(_x,_y,_z,t,_u, _w);
  }
} ;



template <class F, class X, class Y, class Z, class U, class W>
inline Fix12356<F,X,Y,Z,U,W> fix12356(F  f, X x, Y y, Z z, U u, W w) 
{
  return Fix12356<F,X,Y,Z,U,W>(f,x,y,z,u,w);
}

template <class F, class X, class Y, class Z, class T, class W>
class Fix12346
{ F _f;
  X _x;
  Y _y;
  Z _z;
  T _t;
  W _w;
  public:
  Fix12346(F  f, X x, Y y, Z z, T t, W w) :_f(f),_x(x),_y(y),_z(z),_t(t), _w(w){}
  template <class U >
  double operator()(U u)
  { return _f(_x,_y,_z,_t,u, _w);
  }
} ;



template <class F, class X, class Y, class Z, class T, class W>
inline Fix12346<F,X,Y,Z,T,W> fix12356(F  f, X x, Y y, Z z, T t, W w) 
{
  return Fix12346<F,X,Y,Z,T,W>(f,x,y,z,t,w);
}

template <class F, class X, class Y, class Z, class T, class U>
class Fix12345
{ F _f;
  X _x;
  Y _y;
  Z _z;
  T _t;
  U _u;
  public:
  Fix12345(F  f, X x, Y y, Z z, T t, U u) :_f(f),_x(x),_y(y),_z(z),_t(t), _u(u){}
  template <class W>
  double operator()(W w)
  { return _f(_x,_y,_z,_t,_u, w);
  }
} ;


template <class F, class X, class Y, class Z, class T, class U>
inline Fix12345<F,X,Y,Z,T,U> fix12356(F  f, X x, Y y, Z z, T t, U u) 
{
  return Fix12345<F,X,Y,Z,T,U>(f,x,y,z,t,u);
}


#include <fstream>

template <class F>
void plot(F f, double x1, double x2, const char * filename , int npoints=100, double xunit=1, double yunit=1)
{using namespace std;
  ofstream out(filename);
  for(int i=0;i<=npoints;i++)
    {double x=(x1*(npoints-i)+x2*i)/npoints;
     double y=f(x);
     out<<x/xunit<<'\t'<<y/yunit<<endl;
    } 
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#endif
