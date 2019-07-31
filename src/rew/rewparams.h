#ifndef REWPARAMS_H
#define REWPARAMS_H

#include "../params.h"
#include <iostream>


struct RewParam
{
	double val;
	double def;
	//double initial;
	double twk;
	double errPlus;
	double errMinus;
	const string name;
	const string engine;
	bool active;

	RewParam(string name0="", double def0=0, 
	         double errPlus0=0, double errMinus0=0, string engine0="")
	  	:val(def0), def(def0), twk(0),
	  	 errPlus(errPlus0), errMinus(errMinus0), 
	  	 name(name0), engine(engine0), 
	  	 active(false)
	  	 {
	  	 }
	
	double operator()(){ return val;}

	void setErrors(double plus, double minus){ errPlus=plus; errMinus=minus;}

	void reset(){ val=def; twk=0;}
	void set(double v){ val=def=v; twk=0;}

	RewParam & setTwk(double twk0)
	{
		twk=twk0; 
		val=def*(1+twk*(twk>0?errPlus:errMinus));
		return *this;
	}

};

/// global structure with all variables that can be reweighted 
/// ALL params and variables that need reweighting should be placed 
/// in RewParams structure and accessed directly from global variable
/// rew (e.g. rew.qelNorm.val or rew.qel_cc_vector_mass.val )

struct RewParams
{
	RewParam rEnd; // end of reverse iterator
	
	//	Norm
	RewParam qelNorm;
	RewParam resNorm;
	RewParam disNorm;
	RewParam cohNorm;
	RewParam mecNorm;
	RewParam ccNorm;
	RewParam ncNorm;
	RewParam antyNorm;
	
	RewParam dynNorm0;
	RewParam dynNorm1;
	RewParam dynNorm2;
	RewParam dynNorm3;
	RewParam dynNorm4;
	RewParam dynNorm5;
	RewParam dynNorm6;
	RewParam dynNorm7;
	RewParam dynNorm8;
	RewParam dynNorm9;
	
	// QEL
	RewParam bba07_AEp1;
	RewParam bba07_AEp2;
	RewParam bba07_AEp3;
	RewParam bba07_AEp4;
	RewParam bba07_AEp5;
	RewParam bba07_AEp6;
	RewParam bba07_AEp7;

	RewParam bba07_AMp1;
	RewParam bba07_AMp2;
	RewParam bba07_AMp3;
	RewParam bba07_AMp4;
	RewParam bba07_AMp5;
	RewParam bba07_AMp6;
	RewParam bba07_AMp7;

	RewParam bba07_AEn1;
	RewParam bba07_AEn2;
	RewParam bba07_AEn3;
	RewParam bba07_AEn4;
	RewParam bba07_AEn5;
	RewParam bba07_AEn6;
	RewParam bba07_AEn7;

	RewParam bba07_AMn1;
	RewParam bba07_AMn2;
	RewParam bba07_AMn3;
	RewParam bba07_AMn4;
	RewParam bba07_AMn5;
	RewParam bba07_AMn6;
	RewParam bba07_AMn7;

	RewParam bba07_AAx1;
	RewParam bba07_AAx2;
	RewParam bba07_AAx3;
	RewParam bba07_AAx4;
	RewParam bba07_AAx5;
	RewParam bba07_AAx6;
	RewParam bba07_AAx7;
	
	RewParam zexp_nterms;
	RewParam zexp_q4limit;
	RewParam zexp_tc;
	RewParam zexp_t0;

	RewParam zexp_a0;
	RewParam zexp_a1;
	RewParam zexp_a2;
	RewParam zexp_a3;
	RewParam zexp_a4;
	RewParam zexp_a5;
	RewParam zexp_a6;
	RewParam zexp_a7;
	RewParam zexp_a8;
	RewParam zexp_a9;

	RewParam delta_s;
	
	RewParam qel_cc_vector_mass;
	RewParam qel_cc_axial_mass;
	RewParam qel_nc_axial_mass;
	RewParam qel_s_axial_mass;
	RewParam qel_axial_2comp_gamma;
	RewParam qel_axial_2comp_alpha;
	RewParam qel_axial_3comp_theta;
	RewParam qel_axial_3comp_beta;
	
	// RES
	RewParam   pion_axial_mass;  //MaRES
    RewParam   pion_C5A;         //CA5    
    RewParam   SPPBkgScale; 

	
	RewParam End; // end of iferator

	public:
	RewParams():rEnd(),
		
		qelNorm("qelNorm" ,1,0.15,0.10,"rewDyn"),
		resNorm("resNorm" ,1,0,0,"rewDyn"),
		disNorm("disNorm" ,1,0,0,"rewDyn"),
		cohNorm("cohNorm" ,1,0,0,"rewDyn"),
		mecNorm("mecNorm" ,1,0.25,0.25,"rewDyn"),
		ccNorm ("ccNorm"   ,1,0,0,"rewDyn"),
		ncNorm ("ncNorm"   ,1,0,0,"rewDyn"),
		antyNorm("antyNorm",1,0,0,"rewDyn"),

		dynNorm1("dynNorm1",1,0,0,"rewDyn"),
		dynNorm2("dynNorm2",1,0,0,"rewDyn"),
		dynNorm0("dynNorm0",1,0,0,"rewDyn"),
		dynNorm3("dynNorm3",1,0,0,"rewDyn"),
		dynNorm4("dynNorm4",1,0,0,"rewDyn"),
		dynNorm5("dynNorm5",1,0,0,"rewDyn"),
		dynNorm6("dynNorm6",1,0,0,"rewDyn"),
		dynNorm7("dynNorm7",1,0,0,"rewDyn"),
		dynNorm8("dynNorm8",1,0,0,"rewDyn"),
		dynNorm9("dynNorm9",1,0,0,"rewDyn"),


		bba07_AEp1("bba07_AEp1", 1.0, 0.1, 0.1, "rewQEL"),
		bba07_AEp2("bba07_AEp2", 0.9927, 0.1, 0.1, "rewQEL"),
		bba07_AEp3("bba07_AEp3", 0.9898, 0.1, 0.1, "rewQEL"),
		bba07_AEp4("bba07_AEp4", 0.9975, 0.1, 0.1, "rewQEL"),
		bba07_AEp5("bba07_AEp5", 0.9812, 0.1, 0.1, "rewQEL"),
		bba07_AEp6("bba07_AEp6", 0.9340, 0.1, 0.1, "rewQEL"),
		bba07_AEp7("bba07_AEp7", 1.0, 0.1, 0.1, "rewQEL"),

		bba07_AMp1("bba07_AMp1", 1.0, 0.1, 0.1, "rewQEL"),
		bba07_AMp2("bba07_AMp2", 1.0011, 0.1, 0.1, "rewQEL"),
		bba07_AMp3("bba07_AMp3", 0.9992, 0.1, 0.1, "rewQEL"),
		bba07_AMp4("bba07_AMp4", 0.9974, 0.1, 0.1, "rewQEL"),
		bba07_AMp5("bba07_AMp5", 1.0010, 0.1, 0.1, "rewQEL"),
		bba07_AMp6("bba07_AMp6", 1.0003, 0.1, 0.1, "rewQEL"),
		bba07_AMp7("bba07_AMp7", 1.0, 0.1, 0.1, "rewQEL"),

		bba07_AEn1("bba07_AEn1", 1.0, 0.1, 0.1, "rewQEL"),
		bba07_AEn2("bba07_AEn2", 1.011, 0.1, 0.1, "rewQEL"),
		bba07_AEn3("bba07_AEn3", 1.1392, 0.1, 0.1, "rewQEL"),
		bba07_AEn4("bba07_AEn4", 1.0203, 0.1, 0.1, "rewQEL"),
		bba07_AEn5("bba07_AEn5", 1.1093, 0.1, 0.1, "rewQEL"),
		bba07_AEn6("bba07_AEn6", 1.5429, 0.1, 0.1, "rewQEL"),
		bba07_AEn7("bba07_AEn7", 0.9706, 0.1, 0.1, "rewQEL"),

		bba07_AMn1("bba07_AMn1", 1.0, 0.1, 0.1, "rewQEL"),
		bba07_AMn2("bba07_AMn2", 0.9958, 0.1, 0.1, "rewQEL"),
		bba07_AMn3("bba07_AMn3", 0.9877, 0.1, 0.1, "rewQEL"),
		bba07_AMn4("bba07_AMn4", 1.0193, 0.1, 0.1, "rewQEL"),
		bba07_AMn5("bba07_AMn5", 1.0350, 0.1, 0.1, "rewQEL"),
		bba07_AMn6("bba07_AMn6", 0.9164, 0.1, 0.1, "rewQEL"),
		bba07_AMn7("bba07_AMn7", 0.7300, 0.1, 0.1, "rewQEL"),

		bba07_AAx1("bba07_AAx1", 1.0, 0.1, 0.1, "rewQEL"),
		bba07_AAx2("bba07_AAx2", 0.9958, 0.1, 0.1, "rewQEL"),
		bba07_AAx3("bba07_AAx3", 0.9877, 0.1, 0.1, "rewQEL"),
		bba07_AAx4("bba07_AAx4", 1.0193, 0.1, 0.1, "rewQEL"),
		bba07_AAx5("bba07_AAx5", 1.0350, 0.1, 0.1, "rewQEL"),
		bba07_AAx6("bba07_AAx6", 0.9164, 0.1, 0.1, "rewQEL"),
		bba07_AAx7("bba07_AAx7", 0.7300, 0.1, 0.1, "rewQEL"),
		
		zexp_nterms("zexp_nterms", 4, 0, 0, ""), // nonreweightable
		zexp_q4limit("zexp_q4limit", 1, 0, 0, ""), // nonreweightable

		zexp_tc("zexp_tc", 0.1764, 0.1, 0.1, "rewQEL"),
		zexp_t0("zexp_t0", -0.280, 0.1, 0.1, "rewQEL"),

		zexp_a0("zexp_a0", 0.00, 0.1, 0.1, "rewQEL"),
		zexp_a1("zexp_a1", 2.30, 0.1, 0.1, "rewQEL"),
		zexp_a2("zexp_a2", -0.6, 0.1, 0.1, "rewQEL"),
		zexp_a3("zexp_a3", -3.8, 0.1, 0.1, "rewQEL"),
		zexp_a4("zexp_a4", 2.3, 0.1, 0.1, "rewQEL"),
		zexp_a5("zexp_a5", 0.00, 0.1, 0.1, "rewQEL"),
		zexp_a6("zexp_a6", 0.00, 0.1, 0.1, "rewQEL"),
		zexp_a7("zexp_a7", 0.00, 0.1, 0.1, "rewQEL"),
		zexp_a8("zexp_a8", 0.00, 0.1, 0.1, "rewQEL"),
		zexp_a9("zexp_a9", 0.00, 0.1, 0.1, "rewQEL"),

		delta_s("delta_s", -0.15, 0.1, 0.1, "rewQEL"),
		
		qel_cc_vector_mass("qel_cc_vector_mass", 840, 0.16, 0.16, "rewQEL"),
		qel_cc_axial_mass("qel_cc_axial_mass", 1200, 0.16, 0.16, "rewQEL"),
		qel_nc_axial_mass("qel_nc_axial_mass", 1350, 0.16, 0.16, "rewQEL"),
		qel_s_axial_mass("qel_s_axial_mass", 1200, 0.16, 0.16, "rewQEL"),
		qel_axial_2comp_gamma("qel_axial_2comp_gamma", 0.15, 0.1, 0.1, "rewQEL"),
		qel_axial_2comp_alpha("qel_axial_2comp_alpha", 2.0, 0.1, 0.1, "rewQEL"),
		qel_axial_3comp_theta("qel_axial_3comp_theta", 0.15, 0.1, 0.1, "rewQEL"),
		qel_axial_3comp_beta("qel_axial_3comp_beta", 2.0, 0.1, 0.1, "rewQEL"),


	    // RES
	    pion_axial_mass("pion_axial_mass", 0.94, 0.1, 0.1, "rewRES"),
        pion_C5A   ("pion_C5A", 0.19, 0.1, 0.1, "rewRES"),
        SPPBkgScale("SPPBkgScale" ,1, 0.26, 0.26,"rewRES"),
	
		End()	
		{
		}

	RewParam * begin(){return &rEnd+1;}

	RewParam * end(){return &End;}

	RewParam & operator()(string name)
	{
		for(RewParam* p=begin(); p!=end(); p++)
		{
			if(p->name==name)
				return *p;
		}
		return End;
	}

	RewParams & reset()
	{
		for(RewParam* p=begin(); p!=end(); p++)
			p->reset();
		return *this;
	}

	void list(ostream& o)
	{
		for(RewParam* p=begin(); p!=end(); p++)
			if(p->engine!="")
				o<<"\t"<<p->name<<"\n";
	}

	void list_vals(ostream& o)
	{
		for(RewParam* p=begin(); p!=end(); p++)
			if(p->engine!="")
				o<<"\t"<<p->name<<"="<<p->val<<"\n";
	}

	RewParams &init(params &p)
	{
		delta_s.set(p.delta_s);		
		qel_cc_vector_mass.set(p.qel_cc_vector_mass);  // TODO: remove?
		qel_cc_axial_mass.set(p.qel_cc_axial_mass);
		qel_nc_axial_mass.set(p.qel_nc_axial_mass);
		qel_s_axial_mass.set(p.qel_s_axial_mass);
		pion_axial_mass.set(p.pion_axial_mass);
		pion_C5A.set(p.pion_C5A);  
	
		return *this;
	}

};

extern RewParams rew;

#endif
