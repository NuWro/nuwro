/*
 * CGeneralTab.cpp
 *
 *  Created on: 19-06-2011
 *      Author: Marcin Boczulak
 */

#include "CGeneralTab.h"
#include <QtGui>




CGeneralTab::CGeneralTab(map<string, string>& params, QWidget *parent) :
	paramsMap(params) {
	numberOfEvents = new QLineEdit();
	numberOfEvents->setText(paramsMap["number_of_events"].c_str());
	connect(numberOfEvents, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( numberOfEventsChanged( const QString& ) ));

	numberOfTestEvents = new QLineEdit();
	numberOfTestEvents->setText(paramsMap["number_of_test_events"].c_str());

	connect(numberOfTestEvents, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( numberOfTestEventsChanged( const QString& ) ));

	randomSeed = new QLineEdit();
	randomSeed->setText(paramsMap["random_seed"].c_str());
	connect(randomSeed, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( randomSeedChanged( const QString& ) ));

	dyn_qel_cc = new QCheckBox();
	dyn_qel_cc->setChecked(QString(paramsMap["dyn_qel_cc"].c_str()).toInt());
	connect(dyn_qel_cc, SIGNAL(toggled(bool)), this, SLOT( setDynQelCc() ));

	dyn_qel_nc = new QCheckBox();
	dyn_qel_nc->setChecked(QString(paramsMap["dyn_qel_nc"].c_str()).toInt());
	connect(dyn_qel_nc, SIGNAL(toggled(bool)), this, SLOT( setDynQelNc ()));

	dyn_res_cc = new QCheckBox();
	dyn_res_cc->setChecked(QString(paramsMap["dyn_res_cc"].c_str()).toInt());
	connect(dyn_res_cc, SIGNAL(toggled(bool)), this, SLOT( setDynResCc() ));

	dyn_res_nc = new QCheckBox();
	dyn_res_nc->setChecked(QString(paramsMap["dyn_res_nc"].c_str()).toInt());
	connect(dyn_res_nc, SIGNAL(toggled(bool)), this, SLOT( setDynResNc ()));

	dyn_dis_cc = new QCheckBox();
	dyn_dis_cc->setChecked(QString(paramsMap["dyn_dis_cc"].c_str()).toInt());
	connect(dyn_dis_cc, SIGNAL(toggled(bool)), this, SLOT( setDynDisCc() ));

	dyn_dis_nc = new QCheckBox();
	dyn_dis_nc->setChecked(QString(paramsMap["dyn_dis_nc"].c_str()).toInt());
	connect(dyn_dis_nc, SIGNAL(toggled(bool)), this, SLOT( setDynDisNc ()));

	dyn_coh_cc = new QCheckBox();
	dyn_coh_cc->setChecked(QString(paramsMap["dyn_coh_cc"].c_str()).toInt());
	connect(dyn_coh_cc, SIGNAL(toggled(bool)), this, SLOT( setDynCohCc() ));

	dyn_coh_nc = new QCheckBox();
	dyn_coh_nc->setChecked(QString(paramsMap["dyn_coh_nc"].c_str()).toInt());
	connect(dyn_coh_nc, SIGNAL(toggled(bool)), this, SLOT( setDynCohNc ()));

	QVBoxLayout *vBox = new QVBoxLayout;


	QLabel *numberOfEventsLabel = new QLabel("number of events");
	//numberOfEventsLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	vBox->addWidget(numberOfEventsLabel);
	vBox->addWidget(numberOfEvents);

	QLabel *numberOfTestEventsLabel = new QLabel("number of test events");
	//numberOfTestEventsLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	vBox->addWidget(numberOfTestEventsLabel);
	vBox->addWidget(numberOfTestEvents);

	QLabel *randomSeedLabel = new QLabel("random seed");
	//randomSeedLabel->setFrameStyle(QFrame::Panel | QFrame::Sunken);
	vBox->addWidget(randomSeedLabel);
	vBox->addWidget(randomSeed);
	// vBox->addStretch(2);
	//vBox->addSpacing(2);

	QGridLayout *gridLayot = new QGridLayout;
	gridLayot->addWidget(new QLabel("CC"), 0, 1);
	gridLayot->addWidget(new QLabel("NC"), 0, 2);
	gridLayot->addWidget(new QLabel("GEL"), 1, 0);
	gridLayot->addWidget(dyn_qel_cc, 1, 1);
	gridLayot->addWidget(dyn_qel_nc, 1, 2);
	gridLayot->addWidget(new QLabel("RES"), 2, 0);
	gridLayot->addWidget(dyn_res_cc, 2, 1);
	gridLayot->addWidget(dyn_res_nc, 2, 2);
	gridLayot->addWidget(new QLabel("DIS"), 3, 0);
	gridLayot->addWidget(dyn_dis_cc, 3, 1);
	gridLayot->addWidget(dyn_dis_nc, 3, 2);
	gridLayot->addWidget(new QLabel("COH"), 4, 0);
	gridLayot->addWidget(dyn_coh_cc, 4, 1);
	gridLayot->addWidget(dyn_coh_nc, 4, 2);

	QVBoxLayout *vBoxDynamic=new QVBoxLayout();
	gridLayot->setAlignment(Qt::AlignTop);
	vBoxDynamic->addLayout(gridLayot);
	QSpacerItem *spacer=new QSpacerItem(100,110);
	vBoxDynamic->addSpacerItem(spacer);
	QGroupBox *dynamicGroup = new QGroupBox(tr("Dynamics:"));
	dynamicGroup->setLayout(vBoxDynamic);

	// layout->addStretch(1);
	// QVBoxLayout *vbox = new QVBoxLayout();
	//vbox->addLayout(layout);
	QHBoxLayout *hBox = new QHBoxLayout;
	vBox->setAlignment(Qt::AlignTop);
	hBox->addWidget(dynamicGroup);
	hBox->addLayout(vBox);
	setLayout(hBox);


}

CGeneralTab::~CGeneralTab() {
	// TODO Auto-generated destructor stub
}

void CGeneralTab::numberOfEventsChanged(const QString& qs) {
	string number_of_events = qs.toUtf8().constData();
	paramsMap["number_of_events"] = number_of_events;
}

void CGeneralTab::numberOfTestEventsChanged(const QString& qs) {
	string number_of_test_events = qs.toUtf8().constData();
	paramsMap["number_of_test_events"] = number_of_test_events;
}

void CGeneralTab::randomSeedChanged(const QString& qs) {
	string random_seed = qs.toUtf8().constData();
	paramsMap["random_seed"] = random_seed;
}

void CGeneralTab::setDynQelCc() {
	if (dyn_qel_cc->checkState()) {
		paramsMap["dyn_qel_cc"] = "1";
	} else {
		paramsMap["dyn_qel_cc"] = "0";
	}

}
void CGeneralTab::setDynQelNc() {
	if (dyn_qel_nc->checkState()) {
		paramsMap["dyn_qel_nc"] = "1";
	} else {
		paramsMap["dyn_qel_nc"] = "0";
	}

}

void CGeneralTab::setDynResCc() {
	if (dyn_res_cc->checkState()) {
		paramsMap["dyn_res_cc"] = "1";
	} else {
		paramsMap["dyn_res_cc"] = "0";
	}

}

void CGeneralTab::setDynResNc() {
	if (dyn_res_nc->checkState()) {
		paramsMap["dyn_res_nc"] = "1";
	} else {
		paramsMap["dyn_res_nc"] = "0";
	}

}

void CGeneralTab::setDynDisCc() {
	if (dyn_dis_cc->checkState()) {
		paramsMap["dyn_dis_cc"] = "1";
	} else {
		paramsMap["dyn_dis_cc"] = "0";
	}

}

void CGeneralTab::setDynDisNc() {
	if (dyn_dis_nc->checkState()) {
		paramsMap["dyn_dis_nc"] = "1";
	} else {
		paramsMap["dyn_dis_nc"] = "0";
	}

}

void CGeneralTab::setDynCohCc() {
	if (dyn_coh_cc->checkState()) {
		paramsMap["dyn_coh_cc"] = "1";
	} else {
		paramsMap["dyn_coh_cc"] = "0";
	}

}

void CGeneralTab::setDynCohNc() {
	if (dyn_coh_nc->checkState()) {
		paramsMap["dyn_coh_nc"] = "1";
	} else {
		paramsMap["dyn_coh_nc"] = "0";
	}

}
