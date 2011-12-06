#ifndef _rotation_
#define _rotation_
#include "vec.h"
#include "vect.h"
#include "quat.h"


class rotation
{ double xx,xy,xz,
         yx,yy,yz,
	 zx,zy,zz;
friend class lorentz;  
public: 
     rotation():xx(1),xy(0),xz(0),
                yx(0),yy(1),yz(0),
		zx(0),zy(0),zz(1){}
     rotation(double a1, double a2, double a3, double b1, double b2, double b3, double c1, double c2, double c3):
                xx(a1),xy(a1),xz(a3),
                yx(b1),yy(b2),yz(b3),
		zx(c1),zy(c2),zz(c3){}
		
     /// rotation in the plane containing a,b by angle(a,b)
     /// a,b MUAT be normalized
     rotation(vec a,vec b)
     {vec c=vecprod(a,b);
      double angle=atan2(c.length(),a*b);
      *this=rotation(c,angle);
     }
     
     /// rotation about axis by angle alfa 
     /// axis MUST be normalized
     rotation(vec axis, double alfa)
     {double s=sin(alfa),c=cos(alfa),c1=1-c;
      axis=axis.dir();
      double x=axis.x,y=axis.y,z=axis.z;
      xx = c+c1*x*x;    xy = c1*x*y-s*z;  xz = c1*x*z+s*y;
      yx = c1*y*x+s*z;  yy = c+c1*y*y;    yz = c1*y*z-s*x;
      zx = c1*z*x-s*y;  zy = c1*z*y+s*x;  zz = c+c1*z*z;
     }
     
     /// rotation defined by the quaternion h
     rotation (quat h)
     {double x=h.x, y=h.y, z=h.z, w=h.w;
      xx = 1-2*y*y-2*z*z; xy = 2*(x*y-w*z);   xz = 2*(z*x+w*y);
      yx = 2*(x*y+w*z);   yy = 1-2*x*x-2*z*z; yz = 2*(y*z-w*x);
      zx = 2*(z*x-w*y);   zy = 2*(y*z+w*x);   zz = 1-2*x*x-2*y*y;
     }
     
     
     
friend vec operator*(rotation r,vec a) 
{  return vec(r.xx*a.x+r.xy*a.y+r.xz*a.z,
              r.yx*a.x+r.yy*a.y+r.yz*a.z,
              r.zx*a.x+r.zy*a.y+r.zz*a.z);
}

friend vect operator*(rotation r,vect a) 
{  return vect(a.t,
              r.xx*a.x+r.xy*a.y+r.xz*a.z,
              r.yx*a.x+r.yy*a.y+r.yz*a.z,
              r.zx*a.x+r.zy*a.y+r.zz*a.z);
}

/// rotation 
friend rotation operator*(rotation a, rotation b)
{  return rotation(a.xx*b.xx+a.xy*b.yx-a.xz*b.zx,
                   a.xx*b.xy+a.xy*b.yy-a.xz*b.zy,
                   a.xx*b.xz+a.xy*b.yz-a.xz*b.zz,
                   a.yx*b.xx+a.yy*b.yx-a.yz*b.zx,
                   a.yx*b.xy+a.yy*b.yy-a.yz*b.zy,
                   a.yx*b.xz+a.yy*b.yz-a.yz*b.zz,
                   a.zx*b.xx+a.zy*b.yx-a.zz*b.zx,
                   a.zx*b.xy+a.zy*b.yy-a.zz*b.zy,
                   a.zx*b.xz+a.zy*b.yz-a.zz*b.zz);

}

rotation inv()
{return rotation(xx,yx,zx,xy,yy,zy,xz,yx,zz);
}


};

#endif
