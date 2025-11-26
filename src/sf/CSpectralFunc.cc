#include "CSpectralFunc.h"

#include "CgridSpectralFunc.h"
#include "CGaussianSF.h"
using namespace std;

class sftable
{
 CSpectralFunc* t[targEnd][mdEnd][isoEnd];
 
 public: 
 CSpectralFunc* &f(const TargetElement i_target,
					    const MomDistribs i_momDistrib,
					    const IsospinOfSF i_isospin)
 {// cout<<"sf("<<i_target<<","<<i_momDistrib<<","<<i_isospin<<")..."<<endl; 	 
	 return t[i_target][i_momDistrib][i_isospin];
 }
 
 sftable ()
 {
 	 //cout<<"sf init"<<endl;
   for( int k1=0;k1<targEnd;k1++)
	for( int k2=0;k2<mdEnd;k2++)
    for( int k3=0;k3<isoEnd;k3++)
      t[k1][k2][k3]=0;	 
  
  }
 ~sftable (){
    for( int k1=0;k1<targEnd;k1++)
	for( int k2=0;k2<mdEnd;k2++)
    for( int k3=0;k3<isoEnd;k3++)
      {if( t[k1][k2][k3])
  	   //cout<<"sf destroing "<<k1<<","<<k2<<","<<k3<<"..."<<endl; 	 
      delete t[k1][k2][k3];	 
      }
	 //cout<<"sf destroyed"<<endl;
	}
};


using namespace std;
CSpectralFunc* createSF(const TargetElement i_target,
					    const MomDistribs i_momDistrib,
					    const IsospinOfSF i_isospin)
{
    static sftable a;
    CSpectralFunc* &f1=a.f(i_target,i_momDistrib,i_isospin);
	//cout<< "f1="<<f1<<endl;
//	cout<<" creating "<<i_target<<','<<i_momDistrib<<','<<i_isospin<<endl<<endl;
	if(f1) return f1;
//	cout<<"dalej creating "<<i_target<<','<<i_momDistrib<<','<<i_isospin<<endl<<endl;
	
	switch (i_target)
	{
			
		case targC_Ben:
			return f1=
			        new CgridSpectralFunc(
	                   "sf/pke12_tot.grid",
		                createDistrib(md_C12_Ben, proton),
                        1.0/carbon.Z,
			            carbon.PBlock,
			            carbon.Beta 
						);
		case targO_Ben:
			return f1=
			       new CgridSpectralFunc(
	                   "sf/pke16.grid",
		                createDistrib(md_O16_Ben, proton),
                        1.0/oxygen.A,
			            oxygen.PBlock,
			            oxygen.Beta 
						);
						
		case targAr_Ben:
			return f1=
			        new CgridSpectralFunc(
	                    (i_isospin == proton 
					      ? "sf/gsf_Ar40P.grid" 
						  : "sf/gsf_Ar40N.grid"
						),
		                createDistrib(md_Ca40_GCo, i_isospin),
                        1.0,
			            argon.PBlock,
			            argon.Beta
						);
						
		case targFe_Ben:
			return f1=
			        new CgridSpectralFunc(
	                   "sf/pke56_tot.grid",
		                createDistrib(md_Fe56_Ben, proton),
                        1.0771/iron.Z,
			            iron.PBlock,
			            iron.Beta
						);
			
		case targO_GSF:
			return f1=
			        new CGaussianSF(i_momDistrib, i_isospin,
			            oxygen,
						0,5/fermi,0,5/fermi,OxygenP,OxygenN);

		case targCa_GSF:
			return f1=
			        new CGaussianSF(i_momDistrib, i_isospin,
			            calcium,
						0,5/fermi,0,5/fermi,CalciumP,CalciumN);

		case targAr_GSF:
			return f1=
			         new CGaussianSF(i_momDistrib, i_isospin,
			            argon,
						0,5/fermi,0,5/fermi,ArgonP,ArgonN);

		default:
			return 0;

	}
}
