#ifndef _geomy_h_
#define _geomy_h_

#include <iostream>
#include <string>			//geomy constructor
#include <vector>			//
#include <algorithm>		//sort, unique
#include <fstream>			//isotops.txt
#include <sstream>			// ^
#include <stdexcept>		//This header defines a set of standard exceptions that both the library and programs can use to report common errors.
							//function setN(material), vector out of range
#include <cmath> 			//for std::abs in d_eqls (comparing doubles)

#include <TROOT.h>			//
#include <TSystem.h>		//
#include <TFile.h>			//loading geometry file
#include <TGeometry.h>		//geo
#include <TGeoManager.h>	
#include <TGeoBBox.h>		//detector bounding box
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoElement.h>
#include "vec.h"			//
#include "dirs.h"

struct material
{
	double A, Z, N, w_density;
	vec r;
//	double e_bin, p_fermi;
	material(double=0, double=0, double=0);
};


inline material::material(double _a, double _z, double _d)
: A(_a), Z(_z), w_density(_d), N(_a - _z)
{}



inline int d_round(double x) 
{
  return (int)(x + 0.5);
}


class geomy
{
	vec orig, dxyz;
	double max_density;
	double max_length;

	TGeometry* geo;
    TGeoVolume* top;
    TGeoShape* sha1;
	TGeoBBox* box1;
	TGeoNode* node1;
	TGeoMaterial* mat1;
	TGeoMixture* mix1;
	TGeoElement* ele1;



public:
    geomy(std::string _filename, std::string _geoname, std::string volume="", vec _d = vec(0,0,0), vec _o = vec(0,0,0) )
    {		
//		mtrand.SetSeed(0);  // If seed is 0 (default value) a TUUID is generated and used to fill/*
							// the first 8 integers of the seed array.
							// In this case the seed is guaranteed to be unique in space and time.
        top=NULL;
        if(!ifstream(_filename.c_str()))
        	_filename=get_data_dir()+_filename;
        TFile f(_filename.c_str());
        geo = (TGeometry*)f.Get(_geoname.c_str());
        if(!geo)
           throw string("Error reading geometry '")+_geoname+"' from from file: '"+_filename+"'.";
        
        if(volume.length()>0)
           top=gGeoManager->FindVolumeFast(volume.c_str(),false);
        if(top==NULL)   
		top = gGeoManager->GetTopVolume();
		if(top==NULL)
		  throw(" Error making detector");
		sha1 = top->GetShape();
		if(sha1)
		{
				vec a,b;
				sha1->GetAxisRange(1,a.x,b.x);
				sha1->GetAxisRange(2,a.y,b.y);
				sha1->GetAxisRange(3,a.z,b.z);
				orig=(a+b)/2;
				dxyz=(b-a)/2;
				cout<< "Top Volume orig= "<<orig <<";  dxyz ="<<dxyz<<";"<<endl;				
		}
		else 
		  {orig=vec(0,0,0);
		   dxyz=vec(2000,2000,2000);
	      }
		if(not( _d.x == 0 || _d.y == 0 || _d.z == 0))
		{   // find intersection of box of interest with detector geometry
			vec a=max(orig-dxyz,_o-_d);
			vec b=min(orig+dxyz,_o+_d);
			if(b.x>a.x && b.y>a.y && b.z>a.z)
			  {orig=(a+b)/2;
			   dxyz=(b-a)/2;
		      }  
		     else
		      {
		       cerr<<"Box of interest outside detector. No events will be generated"<<endl;
		       exit(23);
		      }
		}

		std::cout << "\nBOX:"<<" O: " << orig << " D: " << dxyz << "\n";
//		std::cout << "Obj = " << Obj(top,0) << "\n";
//		max_density = MaxDensity(top);
//		max_density = 0;
//		std::cout << "MaxDensity = " << max_density << "\n";


    }


	material getpoint()
	{
		vec r =  orig - dxyz + 2*vec(frandom()*dxyz.x,frandom()*dxyz.y,frandom()*dxyz.z);
		material tam=getpoint(r);
//		tam.w_length=1;
		return tam;
	}

	material getpoint(vec r)
	{
		node1 = gGeoManager->FindNode(r.x, r.y, r.z);		//	TGeoNode* node1;
		mat1 = node1->GetMedium()->GetMaterial();			//	TGeoMaterial* mat1;
		//if(mat1->GetZ() > 50) mat1->Print();

  	   //std::cout << "Mixture density: " << mat1->GetDensity() << "\n";
		material tam;
		if( mat1->IsMixture() ) 
		{
			mix1 = (TGeoMixture*)mat1;			//	TGeoMixture* mix1;
			int n = mix1->GetNelements();			
	        double *w=mix1->GetWmixt();
	        double s=0;
	        for(int i=0;i<n;i++)
	            s+=w[i];
	        s*=frandom();
			int c = 0;
			double a=0;
			while( (a+=w[c]) < s )
			{
				++c;
			}
			tam.A = mix1->GetAmixt()[c];
			tam.Z = mix1->GetZmixt()[c];
			tam.N = d_round(tam.A) - tam.Z;
		}
		else	//if(ele1 == NULL)
		{
			tam.A = mat1->GetA();
			tam.Z = mat1->GetZ();
			tam.N = d_round(tam.A) - tam.Z;
		}
		tam.w_density = mat1->GetDensity();
		tam.r = r;
		return tam;	
	}


	material getpoint(vec dir, vec start)
	{   static double MaxLen=0;
		double x,dx,y,dy,z,dz;
	    x=y=z=0;
	    dx=dy=dz=1e40;
		if(dir.x!=0)
		{  x=(orig.x-start.x)/dir.x;
		   dx=abs(dxyz.x/dir.x);
	    }
		if(dir.y!=0)
		{ y=(orig.y-start.y)/dir.y;
		  dy=abs(dxyz.y/dir.y);
	     }
		if(dir.z!=0)
		{ z=(orig.z-start.z)/dir.z;
		  dz=abs(dxyz.z/dir.z);
		}
		material tam;

		double ta=max(x-dx,max(y-dy,z-dz));
		double tb=min(x+dx,min(y+dy,z+dz));
		if(tb>ta)
		 {  double len=dir.length()*(tb-ta);
		    if(len>frandom()*MaxLen)
		      {
			    if(len>MaxLen) 
		            MaxLen=len;
			   tam = getpoint(start+ (ta+frandom()*(tb-ta))*dir);
		      }
		    // else tam.density==0 and event will be rejected  
		  }
		return tam;
	}

	~geomy()
	{

	}

private:


	double Obj(TGeoVolume* t, int d=0)
	{   
		double capsum = 0;
		int n = t->GetNdaughters();
		if(n == 0) return t->Capacity();
		else
		{
			for(int i = 0; i < n; i++)
			{   
				TGeoVolume* t1 = t->GetNode(i)->GetVolume();
				double cap = Obj(t1,d+1);
				capsum += cap;
			}
		}
		return capsum;
	}

	double MaxDensity(TGeoVolume* t)
	{   
		double maxden = 0;
		int n = t->GetNdaughters();

		if(n == 0) 
			return t->GetMaterial()->GetDensity();
		else
		{
			for(int i = 0; i < n; i++)
			{ 
				TGeoVolume* t1 = t->GetNode(i)->GetVolume();
				double tmp_den = MaxDensity(t1);
				if( maxden < tmp_den ) 
					maxden = tmp_den;
			}
		}
		return maxden;
	}
};

#endif
