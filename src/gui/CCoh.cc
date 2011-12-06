/*
 * CCoh.cc
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#include "CCoh.h"

CCoh::CCoh(map<string, string>&params, QWidget *parent):
paramsMap(params) {
	//setFrameStyle(QFrame::Panel | QFrame::Sunken);

	coh_mass_correction = new QLineEdit();
	coh_mass_correction->setText(paramsMap["coh_mass_correction"].c_str());
	connect(coh_mass_correction, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( coh_mass_correctionChanged( const QString& ) ));


	QHBoxLayout *hBoxCoh= new QHBoxLayout();
	hBoxCoh->addWidget(new QLabel("coh_mass_correction"));
	hBoxCoh->addWidget(coh_mass_correction);

	QGroupBox *cohGroup = new QGroupBox("Coh");
	cohGroup->setLayout(hBoxCoh);

	QVBoxLayout *vBox=new QVBoxLayout();
	vBox->addWidget(cohGroup);
	setLayout(vBox);

}

CCoh::~CCoh() {
	// TODO Auto-generated destructor stub
}
void CCoh::coh_mass_correctionChanged(const QString&qs)
{
	paramsMap["coh_mass_correction"] = qs.toUtf8().constData();
}
