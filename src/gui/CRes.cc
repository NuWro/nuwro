/*
 * CRes.cc
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#include "CRes.h"

CRes::CRes(map<string, string>& params, QWidget *parent):
paramsMap(params){
	//setFrameStyle(QFrame::Panel | QFrame::Sunken);

	delta_FF_set = new QLineEdit();
	delta_FF_set->setText(paramsMap["delta_FF_set"].c_str());
	connect(delta_FF_set, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( delta_FF_setChanged( const QString& ) ));

	pion_axial_mass = new QLineEdit();
	pion_axial_mass->setText(paramsMap["pion_axial_mass"].c_str());
	connect(pion_axial_mass, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( pion_axial_massChanged( const QString& ) ));

	pion_C5A = new QLineEdit();
	pion_C5A->setText(paramsMap["pion_C5A"].c_str());
	connect(pion_C5A, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( pion_C5AChanged( const QString& ) ));

	spp_precision = new QLineEdit();
	spp_precision->setText(paramsMap["spp_precision"].c_str());
	connect(spp_precision, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( spp_precisionChanged( const QString& ) ));


	res_dis_cut = new QLineEdit();
	res_dis_cut->setText(paramsMap["res_dis_cut"].c_str());
	connect(res_dis_cut, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( res_dis_cutChanged( const QString& ) ));



	QGridLayout *gridLayot = new QGridLayout();
	gridLayot->addWidget(new QLabel("delta_FF_set"),0,0);
	gridLayot->addWidget(delta_FF_set,0,1);
	gridLayot->addWidget(new QLabel("pion_axial_mass"),0,2);
	gridLayot->addWidget(pion_axial_mass,0,3);
	gridLayot->addWidget(new QLabel("pion_C5A"),1,0);
	gridLayot->addWidget(pion_C5A,1,1);
	QGroupBox *resGroup = new QGroupBox("Res");
	resGroup->setLayout(gridLayot);

	QHBoxLayout *hBoxResDis= new QHBoxLayout();
	hBoxResDis->addWidget(new QLabel("spp_precision"));
	hBoxResDis->addWidget(spp_precision);
	hBoxResDis->addWidget(new QLabel("res_dis_cut"));
	hBoxResDis->addWidget(res_dis_cut);
	QGroupBox *disResGroup = new QGroupBox("Res-dis");
	disResGroup->setLayout(hBoxResDis);

	QVBoxLayout *vBox=new QVBoxLayout();
	vBox->addWidget(resGroup);
	vBox->addWidget(disResGroup);
	setLayout(vBox);

}

CRes::~CRes() {
	// TODO Auto-generated destructor stub
}

void CRes::delta_FF_setChanged(const QString&qs)
{
	paramsMap["delta_FF_set"] = qs.toUtf8().constData();
}
void CRes::pion_axial_massChanged(const QString&qs)
{
	paramsMap["pion_axial_mass"] = qs.toUtf8().constData();
}
void CRes::pion_C5AChanged(const QString&qs)
{
	paramsMap["pion_C5A"] = qs.toUtf8().constData();
}
void CRes::spp_precisionChanged(const QString&qs)
{
	paramsMap["spp_precision"] = qs.toUtf8().constData();
}
void CRes::res_dis_cutChanged(const QString&qs)
{
	paramsMap["res_dis_cut"] = qs.toUtf8().constData();
}
