/*
 * CBeamTab.cpp
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#include "CBeamTab.h"

CBeamTab::CBeamTab(map<string, string> &params, QWidget *parent)
:paramsMap(params){

	beamType = new QComboBox();
	beamType->addItem("0");
	beamType->addItem("1");
	beamType->addItem("2");
	beamType->addItem("3");
	//beamType->setItemText(QString(paramsMap["beam_type"].c_str()).toInt(),paramsMap["beam_type"].c_str());
	//beamType->setText(paramsMap["beam_type"].c_str());
	connect(beamType, SIGNAL( activated( const QString& ) ), this,
			SLOT( beamTypeChanged( const QString& ) ));

//beam type 0 ----------------------------------

	beamParticle = new QLineEdit();
	beamParticle->setText(paramsMap["beam_particle"].c_str());
	connect(beamParticle, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( beamParticleChanged( const QString& ) ));

	beamDirection = new QLineEdit();
	beamDirection->setText(paramsMap["beam_direction"].c_str());
	connect(beamDirection, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( beamDirectionChanged( const QString& ) ));

	beamEnergy = new QLineEdit();
	beamEnergy->setText(paramsMap["beam_energy"].c_str());
	connect(beamEnergy, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( beamEnergyChanged( const QString& ) ));
//end beam type 0---------------------

//start beam type 1---------------------------------------------
	beamContent = new QLineEdit();
	beamContent->setText(paramsMap["beam_content"].c_str());
	connect( beamContent, SIGNAL(textChanged( const QString& ) ),
					                 this, SLOT( beamContentChanged( const QString& ) ) );
//end beam type 1---------------------------------------------

//start beam type 1 2 3---------------------------------------------
	beamFileName = new QLineEdit();
	connect( beamFileName, SIGNAL(textChanged( const QString& ) ),
						                 this, SLOT( beamFileNameChanged( const QString& ) ) );
	QPushButton *LoadFile=new QPushButton("...");
	LoadFile->setFixedWidth(25);
	connect(LoadFile, SIGNAL(clicked()), this, SLOT(loadFileName()));
//end beam type 1 2 3---------------------------------------------


	QVBoxLayout *vBox = new QVBoxLayout;
	QHBoxLayout *hBoxBeamType = new QHBoxLayout;
	QHBoxLayout *hBoxLoadFileBeamTypeOne = new QHBoxLayout;
	QHBoxLayout *hBoxLoadFileBeamTypeTwo = new QHBoxLayout;
	QHBoxLayout *hBoxLoadFileBeamTypeThree = new QHBoxLayout;

	QLabel *beamTypeLabel = new QLabel("Beam type:");
	hBoxBeamType->addWidget(beamTypeLabel);
	hBoxBeamType->addWidget(beamType);

	QPoint TopLeft(0,0);
	QPoint BottomRight(30,30);
	QRect SelectionRectangle(TopLeft, BottomRight);
	//QRubberBand outline (QRubberBand::Rectangle);
	hBoxBeamType->setGeometry(SelectionRectangle);

	vBoxBeamTypeZero = new QVBoxLayout;

	QLabel *beamParticleLabel = new QLabel("beam particle");
	vBoxBeamTypeZero->addWidget(beamParticleLabel);
	vBoxBeamTypeZero->addWidget(beamParticle);

	QLabel *beamDirectionLabel = new QLabel("beam direction");
	vBoxBeamTypeZero->addWidget(beamDirectionLabel);
	vBoxBeamTypeZero->addWidget(beamDirection);

	QLabel *beamEnergyLabel = new QLabel("beam energy");
	vBoxBeamTypeZero->addWidget(beamEnergyLabel);
	vBoxBeamTypeZero->addWidget(beamEnergy);
	beamTypeZeroGroup = new QGroupBox("Beam type 0:");
	beamTypeZeroGroup->setLayout(vBoxBeamTypeZero);


	vBoxBeamTypeOne = new QVBoxLayout;
	QLabel *beamFileNameLabel= new QLabel("Beam File:");
	hBoxLoadFileBeamTypeOne->addWidget(beamFileName);
	hBoxLoadFileBeamTypeOne->addWidget(LoadFile);
	vBoxBeamTypeOne->addWidget(beamFileNameLabel);
	vBoxBeamTypeOne->addLayout(hBoxLoadFileBeamTypeOne);
	beamTypeOneGroup = new QGroupBox("Beam type 1,2,3:");
	beamTypeOneGroup->setLayout(vBoxBeamTypeOne);

//	vBoxBeamTypeTwo = new QVBoxLayout;
//	QLabel *beamFileNameLabel= new QLabel("Beam File:");
//	hBoxLoadFileBeamTypeTwo->addWidget(beamFileName);
//	hBoxLoadFileBeamTypeTwo->addWidget(LoadFile);
//	vBoxBeamTypeTwo->addWidget(beamFileNameLabel);
//	vBoxBeamTypeTwo->addLayout(hBoxLoadFileBeamTypeTwo);
//	beamTypeTwoGroup = new QGroupBox("Beam type 2:");
//	beamTypeTwoGroup->setLayout(vBoxBeamTypeTwo);

//	vBoxBeamTypeThree = new QVBoxLayout;
//	QLabel *beamFileNameLabel= new QLabel("Beam File:");
//	hBoxLoadFileBeamTypeThree->addWidget(beamFileName);
//	hBoxLoadFileBeamTypeThree->addWidget(LoadFile);
//	vBoxBeamTypeThree->addWidget(beamFileNameLabel);
//	vBoxBeamTypeThree->addLayout(hBoxLoadFileBeamTypeThree);
//	beamTypeThreeGroup = new QGroupBox("Beam type 3:");
//	beamTypeThreeGroup->setLayout(vBoxBeamTypeThree);


	//vBox->addLayout(vBoxBeamTypeZero);
    hBoxBeamType->setAlignment(Qt::AlignTop);
	vBox->addLayout(hBoxBeamType);
	beamTypeZeroGroup->setAlignment(Qt::AlignTop);
	vBox->addWidget(beamTypeZeroGroup);
	beamTypeOneGroup->setAlignment(Qt::AlignTop);
	vBox->addWidget(beamTypeOneGroup);
	vBox->setAlignment(Qt::AlignTop);
	beamTypeOneGroup->hide();
//	vBox->addWidget(beamTypeTwoGroup);
//	vBox->addWidget(beamTypeThreeGroup);
//	vSpacer=new QSpacerItem(0, 100);
//	vBox->addSpacerItem(vSpacer);
	setLayout(vBox);


}


CBeamTab::~CBeamTab() {
	// TODO Auto-generated destructor stub
}

void CBeamTab::beamTypeChanged(const QString& qs) {
	string beam_type = qs.toUtf8().constData();
	if(beam_type == "0"){
		beamTypeZeroGroup->show();
		paramsMap["beam_type"] = "0";
		beamTypeOneGroup->hide();
//		beamTypeTwoGroup->hide();
//		beamTypeThreeGroup->hide();

//		vBoxBeamTypeZero->setEnabled(false);
	}
	else if(beam_type == "1"){
		beamTypeOneGroup->show();
		paramsMap["beam_type"] = "1";
		beamTypeZeroGroup->hide();
//		beamTypeTwoGroup->hide();
//		beamTypeThreeGroup->hide();
//		vBoxBeamTypeZero->setEnabled(true);
	}
	else if(beam_type == "2"){
		beamTypeOneGroup->show();		//jak bedzie bardziej dokladnie to porawic na two
		paramsMap["beam_type"] = "2";
		beamTypeZeroGroup->hide();
//		beamTypeOneGroup->hide();
//		beamTypeThreeGroup->hide();
//		vBoxBeamTypeZero->setEnabled(true);
	}
	else if(beam_type == "3"){
		beamTypeOneGroup->show();		//poprawic na three
		paramsMap["beam_type"] = "3";
		beamTypeZeroGroup->hide();
//		beamTypeOneGroup->hide();
//		beamTypeTwoGroup->hide();
//		vBoxBeamTypeZero->setEnabled(true);
	}
}

void CBeamTab::beamParticleChanged(const QString& qs) {
	paramsMap["beam_particle"] = qs.toUtf8().constData();
}

void CBeamTab::beamDirectionChanged(const QString& qs) {
	paramsMap["beam_direction"] = qs.toUtf8().constData();
}

void CBeamTab::beamEnergyChanged(const QString& qs) {
	paramsMap["beam_energy"] = qs.toUtf8().constData();
}

void CBeamTab::beamContentChanged(const QString& qs)
{
	paramsMap["beam_content"] = qs.toUtf8().constData();
}

void CBeamTab::loadFileName()
{
     QString fileName = QFileDialog::getOpenFileName(this,
         tr("Open Beam File"), "",
         tr("Beam File (*.txt);;All Files (*)"));
     if (fileName.isEmpty())
         return;
     else {

         QFile file(fileName);

         if (!file.open(QIODevice::ReadOnly)) {
             QMessageBox::information(this, tr("Unable to open file"),
                 file.errorString());
             return;
         }
     }
     paramsMap["beam_content"] =("@" + fileName).toUtf8().constData();
     beamFileName->setText(fileName.toUtf8().constData());
}

void CBeamTab::beamFileNameChanged(const QString& qs)
{
	paramsMap["beam_content"] =("@" + qs).toUtf8().constData();
}
