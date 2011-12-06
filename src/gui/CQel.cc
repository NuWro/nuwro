/*
 * CGel.cc
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#include "CQel.h"

CQel::CQel(map<string, string>&params, QWidget *parent):
paramsMap(params){
	//setFrameStyle(QFrame::Panel | QFrame::Sunken);

//	QSpinBox *qel_kinematics;// 0-4

	qel_new = new QSpinBox();
	qel_new->setValue(QString(paramsMap["qel_new"].c_str()).toInt());
	qel_new->setMinimum(0);
	qel_new->setMaximum(1);
//	qel_new->setButtonSymbols(PlusMinus);
	connect(qel_new, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( qel_newChanged( const QString& ) ));

	qel_cc_axial_mass = new QLineEdit();
	qel_cc_axial_mass->setText(paramsMap["qel_cc_axial_mass"].c_str());
	connect(qel_cc_axial_mass, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( qel_cc_axial_massChanged( const QString& ) ));

	qel_nc_axial_mass = new QLineEdit();
	qel_nc_axial_mass->setText(paramsMap["qel_nc_axial_mass"].c_str());
	connect(qel_nc_axial_mass, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( qel_nc_axial_massChanged( const QString& ) ));

	qel_cc_vector_ff_set = new QSpinBox();
	qel_cc_vector_ff_set->setValue(QString(paramsMap["qel_cc_vector_ff_set"].c_str()).toInt());
	qel_cc_vector_ff_set->setMaximum(2);
	qel_cc_vector_ff_set->setMinimum(1);
	connect(qel_cc_vector_ff_set, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( qel_cc_vector_ff_setChanged( const QString& ) ));

	qel_cc_axial_ff_set = new QSpinBox();
	qel_cc_axial_ff_set->setValue(QString(paramsMap["qel_cc_axial_ff_set"].c_str()).toInt());
	qel_cc_axial_ff_set->setMinimum(1);
	qel_cc_axial_ff_set->setMaximum(4);
	connect(qel_cc_axial_ff_set, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( qel_cc_axial_ff_setChanged( const QString& ) ));

	qel_strange = new QSpinBox();
	qel_strange->setValue(QString(paramsMap["qel_strange"].c_str()).toInt());
	qel_strange->setMaximum(1);
	qel_strange->setMinimum(0);
	connect(qel_strange, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( qel_strangeChanged( const QString& ) ));

	delta_s = new QLineEdit();
	delta_s->setText(paramsMap["delta_s"].c_str());
	connect(delta_s, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( delta_sChanged( const QString& ) ));

	qel_s_axial_mass = new QLineEdit();
	qel_s_axial_mass->setText(paramsMap["qel_s_axial_mass"].c_str());
	connect(qel_s_axial_mass, SIGNAL(textChanged( const QString& ) ), this,
			SLOT( qel_s_axial_massChanged( const QString& ) ));

	flux_correction = new QSpinBox();
	flux_correction->setMaximum(1);
	flux_correction->setMinimum(0);
	flux_correction->setValue(QString(paramsMap["flux_correction"].c_str()).toInt());
	connect(flux_correction, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( flux_correctionChanged( const QString& ) ));

	sf_method = new QSpinBox();
	sf_method->setMaximum(2);
	sf_method->setMinimum(0);
	sf_method->setValue(QString(paramsMap["sf_method"].c_str()).toInt());
	connect(sf_method, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( sf_methodChanged( const QString& ) ));

	sf_form_factors = new QSpinBox();
	sf_form_factors->setMaximum(3);
	sf_form_factors->setMinimum(0);
	sf_form_factors->setValue(QString(paramsMap["sf_form_factors"].c_str()).toInt());
	connect(sf_form_factors, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( sf_form_factorsChanged( const QString& ) ));

	cc_smoothing = new QSpinBox();
	cc_smoothing->setMaximum(1);
	cc_smoothing->setMinimum(0);
	cc_smoothing->setValue(QString(paramsMap["cc_smoothing"].c_str()).toInt());
	connect(cc_smoothing, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( cc_smoothingChanged( const QString& ) ));

	QGridLayout *gridLayot = new QGridLayout();
	gridLayot->addWidget(new QLabel("qel_new"),0,0);
	gridLayot->addWidget(qel_new,0,1);
	gridLayot->addWidget(new QLabel("qel_cc_axial_mass"),0,2);
	gridLayot->addWidget(qel_cc_axial_mass,0,3);
	gridLayot->addWidget(new QLabel("qel_nc_axial_mass"),1,0);
	gridLayot->addWidget(qel_nc_axial_mass,1,1);
	gridLayot->addWidget(new QLabel("qel_cc_vector_ff_set"),1,2);
	gridLayot->addWidget(qel_cc_vector_ff_set,1,3);
	gridLayot->addWidget(new QLabel("qel_cc_axial_ff_set"),2,0);
	gridLayot->addWidget(qel_cc_axial_ff_set,2,1);
	gridLayot->addWidget(new QLabel("qel_strange"),2,2);
	gridLayot->addWidget(qel_strange,2,3);
	gridLayot->addWidget(new QLabel("delta_s"),3,0);
	gridLayot->addWidget(delta_s,3,1);
	gridLayot->addWidget(new QLabel("qel_s_axial_mass"),3,2);
	gridLayot->addWidget(qel_s_axial_mass,3,3);
	gridLayot->addWidget(new QLabel("flux_correction"),4,0);
	gridLayot->addWidget(flux_correction,4,1);
	gridLayot->addWidget(new QLabel("sf_method"),4,2);
	gridLayot->addWidget(sf_method,4,3);
	gridLayot->addWidget(new QLabel("sf_form_factors"),5,0);
	gridLayot->addWidget(sf_form_factors,5,1);
	gridLayot->addWidget(new QLabel("cc_smoothing"),5,2);
	gridLayot->addWidget(cc_smoothing,5,3);
	QGroupBox *qelGroup = new QGroupBox("Qel");
	qelGroup->setLayout(gridLayot);
	QVBoxLayout *vBox=new QVBoxLayout();
	vBox->addWidget(qelGroup);
	setLayout(vBox);

}

CQel::~CQel() {
	// TODO Auto-generated destructor stub
}

void CQel::qel_newChanged(const QString& qs) {
	paramsMap["qel_new"] = qs.toUtf8().constData();
}

void CQel::qel_cc_axial_massChanged(const QString& qs) {
	paramsMap["qel_cc_axial_mass"] = qs.toUtf8().constData();
}

void CQel::qel_nc_axial_massChanged(const QString& qs) {
	paramsMap["qel_nc_axial_mass"] = qs.toUtf8().constData();
}

void CQel::qel_cc_vector_ff_setChanged(const QString& qs) {
	paramsMap["qel_cc_vector_ff_set"] = qs.toUtf8().constData();
}
void CQel::qel_cc_axial_ff_setChanged(const QString& qs) {
	paramsMap["qel_cc_axial_ff_set"] = qs.toUtf8().constData();
}

void CQel::qel_strangeChanged(const QString& qs) {
	paramsMap["qel_strange"] = qs.toUtf8().constData();
}
void CQel::delta_sChanged(const QString& qs) {
	paramsMap["delta_s"] = qs.toUtf8().constData();
}

void CQel::qel_s_axial_massChanged(const QString& qs) {
	paramsMap["qel_s_axial_mass"] = qs.toUtf8().constData();
}
void CQel::flux_correctionChanged(const QString& qs) {
	paramsMap["flux_correction"] = qs.toUtf8().constData();
}

void CQel::sf_methodChanged(const QString& qs) {
	paramsMap["sf_method"] = qs.toUtf8().constData();
}
void CQel::sf_form_factorsChanged(const QString& qs) {
	paramsMap["sf_form_factors"] = qs.toUtf8().constData();
}

void CQel::cc_smoothingChanged(const QString& qs) {
	paramsMap["cc_smoothing"] = qs.toUtf8().constData();
}

