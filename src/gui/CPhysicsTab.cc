/*
 * CPhysicsTab.cc
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#include "CPhysicsTab.h"

CPhysicsTab::CPhysicsTab(map<string, string> &params, QWidget *parent)
:paramsMap(params){

	QTabWidget *tab = new QTabWidget();
	tab->addTab( new CQel(paramsMap), QString("qel") );
	tab->addTab( new CRes(paramsMap), QString("res") );
	tab->addTab( new CCoh(paramsMap), QString("coh") );
	tab->addTab( new CKaskada(paramsMap), QString("kaskada") );
	/*CQel *qel = new CQel(paramsMap);
	CRes *res = new CRes(paramsMap);
	CCoh *coh = new CCoh(paramsMap);
	CKaskada *kaskada = new CKaskada(paramsMap);
	//qel->setEnabled(false);
	QVBoxLayout *vBox = new QVBoxLayout();
	vBox->addWidget(qel);
	vBox->addWidget(res);
	vBox->addWidget(coh);
	vBox->addWidget(kaskada);
*/

	QVBoxLayout *vBox = new QVBoxLayout();
	vBox->addWidget(tab);
	setLayout(vBox);

}

CPhysicsTab::~CPhysicsTab() {
	// TODO Auto-generated destructor stub
}
