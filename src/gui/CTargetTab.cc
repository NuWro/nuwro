/*
 * CTargetTab.cpp
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#include "CTargetTab.h"


CTargetTab::CTargetTab(map<string, string> &params, QWidget *parent)
:paramsMap(params)
{

	//todo z paramsMap[target_content] wczytywac kolejne liczby do vektora
	target_content.push_back("1");
	target_content.push_back("0");
	target_content.push_back("2");
	target_content.push_back("0");
	target_content.push_back("0");
	target_content.push_back("0");
	target_content.push_back("8");
	target_content.push_back("8");
	target_content.push_back("1");
	target_content.push_back("34");
	target_content.push_back("220");
	target_content.push_back("1");

	targetType = new QComboBox();
	targetType->addItem("0");
	targetType->addItem("1");
	targetType->addItem("2");
	connect(targetType, SIGNAL( activated( const QString& ) ), this,
		SLOT( targetTypeChanged( const QString& ) ));



	QVBoxLayout *vBox = new QVBoxLayout();
	QHBoxLayout *hBoxTargetType = new QHBoxLayout();

	QLabel *targetTypeLabel = new QLabel("Target type:");
	hBoxTargetType->addWidget(targetTypeLabel);
	hBoxTargetType->addWidget(targetType);
	hBoxTargetType->setAlignment(Qt::AlignTop);
	vBox->addLayout(hBoxTargetType);


	//target type 0
	CElement *targetType0 = new CElement(paramsMap);

    nucleus_model = new QComboBox();
    nucleus_model->addItem("0");
    nucleus_model->addItem("1");
	connect(nucleus_model, SIGNAL( activated( const QString& ) ), this,
		SLOT( nucleus_modelChanged( const QString& ) ));

	QVBoxLayout *vBoxTargetTypeZero = new QVBoxLayout();
	QHBoxLayout *hBoxNucleusModel = new QHBoxLayout();
	QLabel *nucleus_modelLabel = new QLabel("nucleus model:");
	hBoxNucleusModel->addWidget(nucleus_modelLabel);
	hBoxNucleusModel->addWidget(nucleus_model);
	targetTypeZeroGroup = new QGroupBox("target type 0");
	vBoxTargetTypeZero->addWidget(targetType0);
	vBoxTargetTypeZero->addLayout(hBoxNucleusModel);
	targetTypeZeroGroup->setLayout(vBoxTargetTypeZero);
	//end target type 0

	//start target type 1
	CElement *targetTypeOne=new CElement(paramsMap,target_content);

	QVBoxLayout *vBoxTargetTypeOne = new QVBoxLayout();
	vBoxTargetTypeOne->addWidget(targetTypeOne);
	targetTypeOneGroup = new QGroupBox("target type 1");
	targetTypeOneGroup->setLayout(vBoxTargetTypeOne);
	//end target type 1


	//start target type 2
	geo_file = new QLineEdit();
	geo_file->setText(paramsMap["geo_file"].c_str());
	connect( geo_file, SIGNAL(textChanged( const QString& ) ),
						                 this, SLOT( geo_fileChanged( const QString& ) ) );
	QPushButton *LoadFile=new QPushButton("...");
	LoadFile->setFixedWidth(25);
	connect(LoadFile, SIGNAL(clicked()), this, SLOT(loadFileName()));
	QHBoxLayout *hBoxLoadGeoFile = new QHBoxLayout;
	hBoxLoadGeoFile->addWidget(geo_file);
	hBoxLoadGeoFile->addWidget(LoadFile);

	geo_name = new QLineEdit();
	geo_name->setText(paramsMap["geo_name"].c_str());
	connect( geo_name, SIGNAL(textChanged( const QString& ) ),
						                 this, SLOT( geo_nameChanged( const QString& ) ) );
	QLabel * geo_nameLabel=new QLabel("geo name");
	QHBoxLayout *hBoxGeoName = new QHBoxLayout;
	hBoxGeoName->addWidget(geo_nameLabel);
	hBoxGeoName->addWidget(geo_name);

	geo_o = new QLineEdit();
	geo_o->setText(paramsMap["geo_o"].c_str());
	connect( geo_o, SIGNAL(textChanged( const QString& ) ),
						                 this, SLOT( geo_oChanged( const QString& ) ) );
	QLabel * geo_oLabel=new QLabel("geo o");
	QHBoxLayout *hBoxGeoO = new QHBoxLayout;
	hBoxGeoO->addWidget(geo_oLabel);
	hBoxGeoO->addWidget(geo_o);


	geo_d= new QLineEdit();
	geo_d->setText(paramsMap["geo_d"].c_str());
	connect( geo_d, SIGNAL(textChanged( const QString& ) ),
						                 this, SLOT( geo_dChanged( const QString& ) ) );
	QLabel * geo_dLabel=new QLabel("geo d");
	QHBoxLayout *hBoxGeoD = new QHBoxLayout;
	hBoxGeoD->addWidget(geo_dLabel);
	hBoxGeoD->addWidget(geo_d);

	QVBoxLayout *vBoxTargetTypeTwo = new QVBoxLayout;

	vBoxTargetTypeTwo->addLayout(hBoxLoadGeoFile);
	vBoxTargetTypeTwo->addLayout(hBoxGeoName);
	vBoxTargetTypeTwo->addLayout(hBoxGeoO);
	vBoxTargetTypeTwo->addLayout(hBoxGeoD);
	targetTypeTwoGroup= new QGroupBox("target type 2");
	targetTypeTwoGroup->setLayout(vBoxTargetTypeTwo);
	//end target type 2

	vBox->addWidget(targetTypeZeroGroup);
	vBox->addWidget(targetTypeOneGroup);
	vBox->addWidget(targetTypeTwoGroup);
	targetTypeOneGroup->hide();
	targetTypeTwoGroup->hide();
	vBox->setAlignment(Qt::AlignTop);
	setLayout(vBox);

}

CTargetTab::~CTargetTab() {
	// TODO Auto-generated destructor stub
}

void CTargetTab::targetTypeChanged(const QString& qs) {
	string target = qs.toUtf8().constData();
	if(target == "0"){
		targetTypeZeroGroup->show();
//		paramsMap["target_content"]="";
//		paramsMap["target_type"] = "0";
		targetTypeOneGroup->hide();
		targetTypeTwoGroup->hide();
}
	else if(target == "1"){
		targetTypeOneGroup->show();
//		paramsMap["target_type"] = "1";
		targetTypeZeroGroup->hide();
		targetTypeTwoGroup->hide();

	}
	else if(target == "2"){
		targetTypeTwoGroup->show();
//		paramsMap["target_type"] = "2";
		targetTypeZeroGroup->hide();
		targetTypeOneGroup->hide();
	}

}

void CTargetTab::nucleus_modelChanged(const QString& qs) {
	string model = qs.toUtf8().constData();
	if(model == "0"){
		paramsMap["nucleus_model"] = "0";
	}
	else if(model == "1"){
		paramsMap["nucleus_model"] = "1";
	}

}

void CTargetTab::geo_fileChanged(const QString& qs)
{
	paramsMap["geo_file"] =qs.toUtf8().constData();
}

void CTargetTab::loadFileName()
{
     QString fileName = QFileDialog::getOpenFileName(this,
         tr("Open Beam File"), "",
         tr("Beam File (*.root);;All Files (*)"));
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

     paramsMap["geo_file"] =fileName.toUtf8().constData();
     geo_file->setText(fileName.toUtf8().constData());
}

void CTargetTab::geo_nameChanged(const QString& qs)
{
	paramsMap["geo_name"] =qs.toUtf8().constData();
}

void CTargetTab::geo_oChanged(const QString& qs)
{
	paramsMap["geo_o"] =qs.toUtf8().constData();
}

void CTargetTab::geo_dChanged(const QString& qs)
{
	paramsMap["geo_d"] =qs.toUtf8().constData();
}

