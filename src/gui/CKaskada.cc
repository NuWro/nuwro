/*
 * CKaskada.cc
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#include "CKaskada.h"

CKaskada::CKaskada(map<string, string>& params, QWidget *parent):
paramsMap(params){
	kaskada_on = new QCheckBox();
	kaskada_on->setChecked(QString(paramsMap["kaskada_on"].c_str()).toInt());
	connect(kaskada_on, SIGNAL(toggled(bool)), this, SLOT( kaskada_onChanged() ));

	pauli_blocking = new QCheckBox();
	pauli_blocking->setChecked(QString(paramsMap["pauli_blocking"].c_str()).toInt());
	connect(pauli_blocking, SIGNAL(toggled(bool)), this, SLOT( pauli_blockingChanged() ));

	formation_zone = new QLineEdit();
	formation_zone->setText(paramsMap["formation_zone"].c_str());
	connect(formation_zone, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( formation_zoneChanged( const QString& ) ));

	xsec = new QComboBox();
	xsec->addItem("metropolis");
	xsec->addItem("oset");
	//beamType->setItemText(QString(paramsMap["beam_type"].c_str()).toInt(),paramsMap["beam_type"].c_str());
	//beamType->setText(paramsMap["beam_type"].c_str());
	connect(xsec, SIGNAL( activated( const QString& ) ), this,
			SLOT( xsecChanged( const QString& ) ));

	QVBoxLayout *vBox= new QVBoxLayout();
	vBox->addWidget(new QLabel("kaskada_on"));
	vBox->addWidget(kaskada_on);
	vBox->addWidget(new QLabel("pauli_blocking"));
	vBox->addWidget(pauli_blocking);
	vBox->addWidget(new QLabel("formation_zone"));
	vBox->addWidget(formation_zone);
	vBox->addWidget(new QLabel("xsec"));
	vBox->addWidget(xsec);

	setLayout(vBox);




}

CKaskada::~CKaskada() {
	// TODO Auto-generated destructor stub
}



void CKaskada::kaskada_onChanged() {
	if (kaskada_on->checkState()) {
		paramsMap["kaskada_on"] = "1";
	} else {
		paramsMap["kaskada_on"] = "0";
	}

}

void CKaskada::pauli_blockingChanged() {
	if (pauli_blocking->checkState()) {
		paramsMap["pauli_blocking"] = "1";
	} else {
		paramsMap["pauli_blocking"] = "0";
	}

}

void CKaskada::formation_zoneChanged(const QString& qs) {
	paramsMap["formation_zone"] = qs.toUtf8().constData();
}

void CKaskada::xsecChanged(const QString& qs) {
	string beam_type = qs.toUtf8().constData();
	if(beam_type == "metropolis"){
		paramsMap["xsec"] = "0";

	}
	else if(beam_type == "oset"){
		paramsMap["xsec"] = "1";
	}
}

