#ifndef _geomy_h_
#define _geomy_h_

#include <iostream>
#include <string>		//geomy constructor
#include <vector>		//
#include <algorithm>		//sort, unique
#include <fstream>		//isotops.txt
#include <sstream>		// ^
#include <stdexcept>		//This header defines a set of standard exceptions that both the library and programs can use to report common errors.
				//function setN(material), vector out of range
#include <cmath> 		//for std::abs in d_eqls (comparing doubles)

#include <TROOT.h>		//
#include <TSystem.h>		//
#include <TFile.h>		//loading geometry file
#include <TGeometry.h>		//geo
#include <TGeoManager.h>	
#include <TGeoBBox.h>		//_detector bounding box
#include <TGeoMedium.h>
#include <TGeoMaterial.h>
#include <TGeoElement.h>
#include "vec.h"			//
#include "dirs.h"


struct material
{
	double A, Z, N, w_density;
	vec r;
	material(double=0, double=0, double=0);
};


inline material::material(double _a, double _z, double _d)
: A(_a), Z(_z), N(_a - _z), w_density(_d)
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
	double dens,ndens;
	double pots;
	double npots;
	double nuc[2];


	double length_scale=1;
	double density_convert=1;


	TGeometry* geo;
	TGeoVolume* top;
	TGeoShape* sha1;
	TGeoBBox* box1;
	TGeoNode* node1;
	TGeoMaterial* mat1;
	TGeoMixture* mix1;
	TGeoElement* ele1;


//C Thorpe: reads in length units from params file, used in density calculations
public:
	geomy(std::string _filename, std::string _geoname,std::string geom_length_units, double geom_density_convert, std::string volume="", vec _d = vec(0,0,0), vec _o = vec(0,0,0) )
	{		

	//C Thorpe: Factor to convert density to g/cm3 - required value has been added to ND280 geom.txt file
	density_convert=geom_density_convert;
		
	//set length scaling based on units provided
	if(geom_length_units == "mm") length_scale = 0.1;
	else if(geom_length_units == "cm") length_scale = 1;
	else if(geom_length_units == "m") length_scale = 10;
	else
	std::cout << "Unrecognised geometry length units: " << geom_length_units  << std::endl<< "Use either mm, cm or m" << std::endl;

		


		pots=npots=0;


	  	  max_length=0;
	  	  max_density=0;
			dens=ndens=0;
			nuc[0]=nuc[1]=0;
		//		mtrand.SetSeed(0);  // If seed is 0 (default value) a TUUID is generated and used to fill/*
								// the first 8 integers of the seed array.
								// In this case the seed is guaranteed to be unique in space and time.
		top=NULL;
		if(!ifstream(_filename.c_str()))
			_filename=get_data_dir()+_filename;
		TFile f(_filename.c_str());
		TGeoManager::SetVerboseLevel(0);
		geo = (TGeometry*)f.Get(_geoname.c_str()); // TGeoNavigator writes to cerr without reason (how to block it?)
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
		{
			orig=vec(0,0,0);
			dxyz=vec(2000,2000,2000);
		}
		
		if(not( _d.x == 0 || _d.y == 0 || _d.z == 0))
		{   // find intersection of box of interest with detector geometry
			vec a=max(orig-dxyz,_o-_d);
			vec b=min(orig+dxyz,_o+_d);
			if(b.x>a.x && b.y>a.y && b.z>a.z)
			{
				orig=(a+b)/2;
				dxyz=(b-a)/2;
			}  
			else
			{
				cerr<<"Box of interest outside detector. No events will be generated"<<endl;
				exit(23);
			}
		}
		
		std::cout << "\nBOX:"<<" O: " << orig << " D: " << dxyz << "\n";
    }


	material getpoint()
	{
		vec r =  orig - dxyz + 2*vec(frandom()*dxyz.x,frandom()*dxyz.y,frandom()*dxyz.z);
		return getpoint(r);
	}

	material getpoint(vec r)
	{
		node1 = gGeoManager->FindNode(r.x, r.y, r.z);		//	TGeoNode* node1;
		mat1 = node1->GetMedium()->GetMaterial();			//	TGeoMaterial* mat1;

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

			  	
			for(int i=0;i<n;i++){
			  //  std::cout << mix1->GetAmixt()[i] << "   " << mix1->GetZmixt()[i] << "   " << w[i] << std::endl;
			}
			

			tam.A = mix1->GetAmixt()[c];
			tam.Z = mix1->GetZmixt()[c];
			tam.N = d_round(tam.A) - tam.Z;
		}
		else	//if(ele1 == NULL)
		{
		  //  std::cout << "Not mixture" << std::endl;
			tam.A = mat1->GetA();
			tam.Z = mat1->GetZ();
			tam.N = d_round(tam.A) - tam.Z;
		}
		


		//	double CLHEP_g_cm3=6.24151e+18;
			
		//conversion factor to get density in g/cm3	
		//for ND280 density_convert = 6.24151e+18;
		tam.w_density = mat1->GetDensity()/density_convert;

		dens+=tam.w_density;
		ndens++;
		max_density = max(max_density,tam.w_density);
//		std::cout << nuc[0] << std::endl;
		nuc[0]+=tam.N*tam.w_density/(tam.N+tam.Z);
		nuc[1]+=tam.Z*tam.w_density/(tam.N+tam.Z);
		tam.r = r;
		return tam;	
	}

	bool is_hit_by(vec dir, vec start)
	{   
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

		double ta=max(x-dx,max(y-dy,z-dz));
		double tb=min(x+dx,min(y+dy,z+dz));
		return tb>ta;
    }
    
	material getpoint(vec dir, vec start)
	{   
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
		tam.w_density=0;

		double ta=max(x-dx,max(y-dy,z-dz));
		double tb=min(x+dx,min(y+dy,z+dz));
		if(tb>ta)
		 { 
			//C Thorpe: Check geom units, convert to cm
			 
			double len=dir.length()*(tb-ta)*length_scale; //length scale takes geom length units and converts to cm
	
		    if(len>frandom()*max_length)
		      {
			    if(len>max_length) 
		            max_length=len;
			    tam = getpoint(start+ (ta+frandom()*(tb-ta))*dir);
			    
		      }
		 }
		double mol=6.02214129e23 ;

		/* density is in g/cm3 
		*  length is in cm
		*  mol = number of nucleons per gram
		*  pots is in nucleons/cm2
		*/ 
		//CT: this density is wrong - off by many OM
		// max_length is in cm, tam.w_density is in g/cm3
		
		pots+=tam.w_density*max_length*mol;
		npots++;  
		

		//	std::cout << npots << std::endl;
		//~ cout<<tam.w_density/CLHEP_g_cm3<<'\t';
		//~ cout<<tam.A<<'\t';
		//~ cout<<tam.Z<<'\t';
		//~ 
		//~ cout<<max_density/CLHEP_g_cm3<<'\t';
		//~ cout<<max_length/10<<'\t';
		//~ cout<<nucleons_per_cm2()<<endl;
		return tam;
	}
	
	double max_len()
	{
		return max_length;
	}

	double max_dens()
	{
		return max_density;
	}
	
	double nucleons_per_cm2()
	{
		return npots>0 ? pots/npots : 0;
	}
	
	double density()
	{
		return dens/ndens*gram/cm3; 
	}
	
	double vol_mass()
	{
	 	//C Thorpe: Lengths should be in cm and density in g/cm3 for all configurations 
	 
	return 8*dxyz.x*dxyz.y*dxyz.z*density()*length_scale*length_scale*length_scale*cm3; 

	}

	double frac_proton()
	{
		return nuc[1]/(nuc[0]+nuc[1]);
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
