/*
 * CPhysicsTab.h
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#ifndef CPHYSICSTAB_H_
#define CPHYSICSTAB_H_

#include <QtGui>
#include <map>

#include "CQel.h"
#include "CRes.h"
#include "CCoh.h"
#include "CKaskada.h"

using namespace std;


class CPhysicsTab : public QWidget{
Q_OBJECT
public:
	CPhysicsTab(map<string, string> &, QWidget *parent = 0);
	virtual ~CPhysicsTab();
private:
	map<string, string> &paramsMap;
private slots:
};

#endif /* CPHYSICSTAB_H_ */
