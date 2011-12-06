/*
 * CMainWindow.cc
 *
 *  Created on: 14-06-2011
 *      Author: boczus
 */
#include "CMainWindow.h"
#include "params.h"
using namespace std;

CMainWindow::CMainWindow(){


		setDefaultParams(); //todo zamiana na classe wczytujaca z pliku poczatkowe parametry

		resize(400, 400);
    	setWindowTitle(
            QApplication::translate("QtNuwro", "QtNuwro"));


		startButton= new QPushButton("Start");
		quitButton= new QPushButton("Exit");
		quitButton->setGeometry(25, 15, 250, 175);
		startButton->setGeometry(25, 15, 250, 175);
		connect(startButton, SIGNAL(clicked()),this, SLOT(startNuwro()));
		connect(quitButton, SIGNAL(clicked()),qApp, SLOT(quit()));
		widget=new QWidget();

		tab = new QTabWidget();
		tab->addTab( new CGeneralTab(paramsMap), QString("General") );
		tab->addTab( new CBeamTab(paramsMap), QString("Beam") );
		tab->addTab( new CTargetTab(paramsMap), QString("Target") );
		tab->addTab( new CPhysicsTab(paramsMap), QString("Physics") );

		//tab->addTab( new QLabel( QString("Parametry3"), tab, 0 ), QString("Parametry3") );

	    layout = new QVBoxLayout();
	    layout->addWidget(tab);
	    layout->addWidget(startButton);
	    layout->addWidget(quitButton);


	    widget->setLayout(layout);

//osobna clasa na menu
	    QMenuBar *menu= new QMenuBar();
	    menu->addMenu("&File");
	    menu->addMenu("&Edit");
	    setMenuBar(menu);


//	    setCorner(Qt::TopLeftCorner, Qt::LeftDockWidgetArea);
//	    setCorner(Qt::BottomLeftCorner, Qt::LeftDockWidgetArea);
//	    setCorner(Qt::TopRightCorner, Qt::RightDockWidgetArea);
//	    setCorner(Qt::BottomRightCorner, Qt::RightDockWidgetArea);
	    setCentralWidget(widget);



	}

CMainWindow::~CMainWindow()
{

}
void CMainWindow::setDefaultParams()
{
//	map<string,string> mapa;

	
	paramsMap["number_of_events"]="1000";
	paramsMap["number_of_test_events"]="2000";
	paramsMap["random_seed"]="1";
	paramsMap["dyn_qel_cc"]="1";
	paramsMap["dyn_qel_nc"]="0";
	paramsMap["dyn_res_cc"]="0";
	paramsMap["dyn_res_nc"]="0";
	paramsMap["dyn_dis_cc"]="0";
	paramsMap["dyn_dis_nc"]="0";
	paramsMap["dyn_coh_cc"]="0";
	paramsMap["dyn_coh_nc"]="0";

	paramsMap["beam_type"] = "0";
	paramsMap["beam_particle"] = "14";
	paramsMap["beam_direction"] = "0 0 1";
	paramsMap["beam_energy"] = "350 5550 41 79 119 151 162 169 173 183 150 149 157 145 136 108 92 87 69 70 63 39 45 30 34 37 29 25 21 17 15 12 12 11 10 9 9 9 9 8 7 7 6 6 6 6 6 6 5 5 5 4 4 4 3 3 3";

	paramsMap["beam_content"] = "14 50% 0 3000";

	//target type 0
	paramsMap["nucleus_p"] = "8";
	paramsMap["nucleus_n"] = "8";
	paramsMap["nucleus_E_b"] = "34";
	paramsMap["nucleus_kf"] = "220";
	paramsMap["nucleus_model"] = "0";
	paramsMap["nucleus_target"] = "0";

	//target type 2
	paramsMap["geo_file"] = "target/ND280.root";
	paramsMap["geo_name"] = "ND280Geometry";
	paramsMap["geo_o"] = "0 0 0";
	paramsMap["geo_d"] = "0 0 0";

	//Qel
	paramsMap["qel_new"] = "1";
	paramsMap["qel_cc_axial_mass"] = "1200";
	paramsMap["qel_nc_axial_mass"] = "1200";
	paramsMap["qel_cc_vector_ff_set"] = "2";
	paramsMap["qel_cc_axial_ff_set"] = "1";
	paramsMap["qel_strange"] = "1";
	paramsMap["delta_s"] = "-0.15";
	paramsMap["qel_s_axial_mass"] = "1200";
	paramsMap["flux_correction"] = "0";
	paramsMap["sf_method"] = "0";
	paramsMap["sf_form_factors"] = "1";
	paramsMap["cc_smoothing"] = "0";
//	paramsMap["qel_kinematics"] = "";

	//Res
	paramsMap["delta_FF_set"] = "0.94";
	paramsMap["pion_axial_mass"] = "-0.15";
	paramsMap["pion_C5A"] = "1.19";
//	paramsMap["delta_FF_set"] = "";

	// RES - DIS boundary
	paramsMap["spp_precision"] = "500";
	paramsMap["res_dis_cut"] = "1600";

	// COH
	paramsMap["coh_mass_correction"] = "1";

	//Final state interaction parameters
	//nucleus_model       = 0  //"flatnucleus" i.e. nucleus is a ball
	paramsMap["kaskada_on"] = "1";
//	paramsMap["kaskada_debug"] = "";
	paramsMap["pauli_blocking"] = "1";
	paramsMap["formation_zone"] = "fz";
	paramsMap["first_step"] = "1";
	paramsMap["step"] = "0.2";
	paramsMap["xsec"] = "1";

	#define PARAM(type,name,default_value) {stringstream s; s<<default_value;paramsMap[#name]=s.str();}
    PARAMS_ALL() 
    #undef PARAM

}

void CMainWindow::startNuwro()
{
	argv = new char*[100];
	vector<string> v;
	v.push_back("nuwro");
	map<string,string>::iterator it;
	for(it = paramsMap.begin(); it != paramsMap.end(); ++it)
	{
		if((*it).second[0]!='@')
		{
			v.push_back("-p");
			v.push_back((*it).first + "=" + (*it).second);
		}
		else{
			v.push_back("-p");
			v.push_back((*it).second);
		}

	}
	for (int i = 0;  i < v.size(); ++i)
	{
		argv[i]=new char[v[i].length()+1];
		strcpy(argv[i],v[i].c_str());
	}
	nuwro.main(v.size(),argv);
}


