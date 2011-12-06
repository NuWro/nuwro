/*
 * CElement.cc
 *
 *  Created on: 26-06-2011
 *      Author: boczus
 */

#include "CElement.h"
using namespace std;

CElement::CElement(map<string, string> &params, QWidget *parent)
:paramsMap(params){
	nucleus_p = new QSpinBox();
	nucleus_p->setValue(QString(paramsMap["nucleus_p"].c_str()).toInt());
//	nucleus_p->setButtonSymbols(PlusMinus);
	connect(nucleus_p, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_pChanged( const QString& ) ));

	nucleus_n = new QSpinBox();
	nucleus_n->setValue(QString(paramsMap["nucleus_n"].c_str()).toInt());
	connect(nucleus_n, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_nChanged( const QString& ) ));

	nucleus_E_b = new QLineEdit();
	nucleus_E_b->setText(paramsMap["nucleus_E_b"].c_str());
	connect(nucleus_E_b, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_E_bChanged( const QString& ) ));

	nucleus_kf = new QLineEdit();
	nucleus_kf->setText(paramsMap["nucleus_kf"].c_str());
	connect(nucleus_kf, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_kfChanged( const QString& ) ));


    nucleus_target = new QComboBox();
    nucleus_target->addItem("0");
    nucleus_target->addItem("1");
    nucleus_target->addItem("2");
    nucleus_target->addItem("3");
    nucleus_target->addItem("4");
    nucleus_target->addItem("5");
    nucleus_target->addItem("6");
	connect(nucleus_target, SIGNAL( activated( const QString& ) ), this,
		SLOT( nucleus_targetChanged( const QString& ) ));



	QVBoxLayout *vBox = new QVBoxLayout();
	QHBoxLayout *hBoxTargetType = new QHBoxLayout();

	QGridLayout *gridTargetType = new QGridLayout();
	QLabel *nucleus_pLabel= new QLabel("nucleus p");
	gridTargetType->addWidget(nucleus_pLabel,0,0);
	gridTargetType->addWidget(nucleus_p,0,1);

	QLabel *nucleus_nLabel= new QLabel("nucleus n");
	gridTargetType->addWidget(nucleus_nLabel,1,0);
	gridTargetType->addWidget(nucleus_n,1,1);

	QLabel *nucleus_E_bLabel= new QLabel("nucleus E_b");
	gridTargetType->addWidget(nucleus_E_bLabel,2,0);
	gridTargetType->addWidget(nucleus_E_b,2,1);

	QLabel *nucleus_kfLabel= new QLabel("nucleus kf");
	gridTargetType->addWidget(nucleus_kfLabel,3,0);
	gridTargetType->addWidget(nucleus_kf,3,1);

	QLabel *nucleusTargetLabel = new QLabel("nucleus target:");
	gridTargetType->addWidget(nucleusTargetLabel,5,0);
	gridTargetType->addWidget(nucleus_target,5,1);


	vBox->addLayout(gridTargetType);
	vBox->setAlignment(Qt::AlignTop);
	setLayout(vBox);

}

CElement::CElement (map<string, string> &params,vector<string> content, QWidget *parent)
:paramsMap(params),target_content(content)
 {

	nucleus_p = new QSpinBox();
	nucleus_p->setValue(QString(target_content[0].c_str()).toInt());
	nucleus_p->setButtonSymbols(QSpinBox::QAbstractSpinBox::PlusMinus);
	connect(nucleus_p, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_pTargetContentChanged( const QString& ) ));

	nucleus_n = new QSpinBox();
	nucleus_n->setValue(QString(target_content[1].c_str()).toInt());
	connect(nucleus_n, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_nTargetContentChanged( const QString& ) ));

	count = new QSpinBox();
	count->setValue(QString(target_content[2].c_str()).toInt());
//	nucleus_p->setButtonSymbols(PlusMinus);
	connect(count, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( countTargetContentChanged( const QString& ) ));

	nucleus_E_b = new QLineEdit();
	nucleus_E_b->setText(target_content[3].c_str());
	connect(nucleus_E_b, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_E_bTargetContentChanged( const QString& ) ));

	nucleus_kf = new QLineEdit();
	nucleus_kf->setText(target_content[4].c_str());
	connect(nucleus_kf, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_kfTargetContentChanged( const QString& ) ));


    nucleus_target = new QComboBox();
    nucleus_target->addItem("0");
    nucleus_target->addItem("1");
    nucleus_target->addItem("2");
    nucleus_target->addItem("3");
    nucleus_target->addItem("4");
    nucleus_target->addItem("5");
    nucleus_target->addItem("6");
	connect(nucleus_target, SIGNAL( activated( const QString& ) ), this,
		SLOT( nucleus_targetTargetContentChanged( const QString& ) ));

	nucleus_p2 = new QSpinBox();
	nucleus_p2->setValue(QString(target_content[6].c_str()).toInt());
//	nucleus_p->setButtonSymbols(PlusMinus);
	connect(nucleus_p2, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_p2TargetContentChanged( const QString& ) ));

	nucleus_n2 = new QSpinBox();
	nucleus_n2->setValue(QString(target_content[7].c_str()).toInt());
	connect(nucleus_n2, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( nucleus_n2TargetContentChanged( const QString& ) ));

	count2 = new QSpinBox();
	count2->setValue(QString(target_content[8].c_str()).toInt());
//	nucleus_p->setButtonSymbols(PlusMinus);
	connect(count2, SIGNAL( valueChanged( const QString& ) ), this,
			SLOT( count2TargetContentChanged( const QString& ) ));

	nucleus_E_b2 = new QLineEdit();
	nucleus_E_b2->setText(target_content[9].c_str());
	connect(nucleus_E_b2, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_E_b2TargetContentChanged( const QString& ) ));

	nucleus_kf2 = new QLineEdit();
	nucleus_kf2->setText(target_content[10].c_str());
	connect(nucleus_kf2, SIGNAL( textChanged( const QString& ) ), this,
			SLOT( nucleus_kf2TargetContentChanged( const QString& ) ));


    nucleus_target2 = new QComboBox();
    nucleus_target2->addItem("0");
    nucleus_target2->addItem("1");
    nucleus_target2->addItem("2");
    nucleus_target2->addItem("3");
    nucleus_target2->addItem("4");
    nucleus_target2->addItem("5");
    nucleus_target2->addItem("6");
	connect(nucleus_target2, SIGNAL( activated( const QString& ) ), this,
		SLOT( nucleus_target2TargetContentChanged( const QString& ) ));



	QVBoxLayout *vBox = new QVBoxLayout();
	QHBoxLayout *hBoxTargetType = new QHBoxLayout();

	QGridLayout *gridTargetType = new QGridLayout();

	QLabel *elementOne = new QLabel("element 1");
	QLabel *elementTwo = new QLabel("element 2");
	gridTargetType->addWidget(elementOne,0,1);
	gridTargetType->addWidget(elementTwo,0,2);
	QLabel *nucleus_pLabel= new QLabel("nucleus p");
	gridTargetType->addWidget(nucleus_pLabel,1,0);
	gridTargetType->addWidget(nucleus_p,1,1);
	gridTargetType->addWidget(nucleus_p2,1,2);
	QLabel *nucleus_nLabel= new QLabel("nucleus n");
	gridTargetType->addWidget(nucleus_nLabel,2,0);
	gridTargetType->addWidget(nucleus_n,2,1);
	gridTargetType->addWidget(nucleus_n2,2,2);
	QLabel *countLabel= new QLabel("count");
	gridTargetType->addWidget(countLabel,3,0);
	gridTargetType->addWidget(count,3,1);
	gridTargetType->addWidget(count2,3,2);
	QLabel *nucleus_E_bLabel= new QLabel("nucleus E_b");
	gridTargetType->addWidget(nucleus_E_bLabel,4,0);
	gridTargetType->addWidget(nucleus_E_b,4,1);
	gridTargetType->addWidget(nucleus_E_b2,4,2);
	QLabel *nucleus_kfLabel= new QLabel("nucleus kf");
	gridTargetType->addWidget(nucleus_kfLabel,5,0);
	gridTargetType->addWidget(nucleus_kf,5,1);
	gridTargetType->addWidget(nucleus_kf2,5,2);
	QLabel *nucleusTargetLabel = new QLabel("nucleus target:");
	gridTargetType->addWidget(nucleusTargetLabel,6,0);
	gridTargetType->addWidget(nucleus_target,6,1);
	gridTargetType->addWidget(nucleus_target2,6,2);


	vBox->addLayout(gridTargetType);
	vBox->setAlignment(Qt::AlignTop);
	setLayout(vBox);


}


CElement::~CElement() {
	// TODO Auto-generated destructor stub
}



void CElement::changeParamTargetContent()
{
    string content="";
    for(int i = 0;i < target_content.size();++i){
        content += target_content[i] + " ";
    }
    paramsMap["target_content"] = content;
}


void CElement::nucleus_pChanged(const QString& qs) {
	paramsMap["nucleus_p"] = qs.toUtf8().constData();
}

void CElement::nucleus_nChanged(const QString& qs) {
	paramsMap["nucleus_n"] = qs.toUtf8().constData();
}

void CElement::nucleus_E_bChanged(const QString& qs) {
	paramsMap["nucleus_E_b"] = qs.toUtf8().constData();
}

void CElement::nucleus_kfChanged(const QString& qs)
{
	paramsMap["nucleus_kf"] = qs.toUtf8().constData();
}


void CElement::nucleus_targetChanged(const QString& qs) {
	string target = qs.toUtf8().constData();
	if(target == "0"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "1"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "2"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "3"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "4"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "5"){
		paramsMap["nucleus_target"] = target;
	}
	else if(target == "6"){
		paramsMap["nucleus_target"] = target;
	}

}


void CElement::nucleus_pTargetContentChanged(const QString& qs) {
	target_content[0] = qs.toUtf8().constData();
	changeParamTargetContent();
}

void CElement::nucleus_nTargetContentChanged(const QString& qs) {
	target_content[1] = qs.toUtf8().constData();
	changeParamTargetContent();

}

void CElement::countTargetContentChanged(const QString& qs) {
	target_content[2] = qs.toUtf8().constData();
	target_content[2] +="x";
	changeParamTargetContent();
}


void CElement::nucleus_E_bTargetContentChanged(const QString& qs) {
	target_content[3] = qs.toUtf8().constData();
	changeParamTargetContent();
}

void CElement::nucleus_kfTargetContentChanged(const QString& qs)
{
	target_content[4] = qs.toUtf8().constData();
	changeParamTargetContent();
}


void CElement::nucleus_targetTargetContentChanged(const QString& qs) {
	string target = qs.toUtf8().constData();
	if(target == "0"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "1"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "2"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "3"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "4"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "5"){
		target_content[5] = target;
		changeParamTargetContent();
	}
	else if(target == "6"){
		target_content[5] = target;
		changeParamTargetContent();
	}

}

void CElement::nucleus_p2TargetContentChanged(const QString& qs) {
	target_content[6] = qs.toUtf8().constData();
	changeParamTargetContent();
}

void CElement::nucleus_n2TargetContentChanged(const QString& qs) {
	target_content[7] = qs.toUtf8().constData();
	changeParamTargetContent();
}

void CElement::count2TargetContentChanged(const QString& qs) {
	target_content[8] = qs.toUtf8().constData();
	target_content[8]+="x";
	changeParamTargetContent();
}


void CElement::nucleus_E_b2TargetContentChanged(const QString& qs) {
	target_content[9] = qs.toUtf8().constData();
	changeParamTargetContent();
}

void CElement::nucleus_kf2TargetContentChanged(const QString& qs)
{
	target_content[10] = qs.toUtf8().constData();
	changeParamTargetContent();
}


void CElement::nucleus_target2TargetContentChanged(const QString& qs) {
	string target = qs.toUtf8().constData();
	if(target == "0"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "1"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "2"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "3"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "4"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "5"){
		target_content[11] = target;
		changeParamTargetContent();
	}
	else if(target == "6"){
		target_content[11] = target;
		changeParamTargetContent();
	}

}
