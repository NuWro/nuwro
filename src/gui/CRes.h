/*
 * CRes.h
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#ifndef CRES_H_
#define CRES_H_

#include <QtGui>
#include <map>
using namespace std;

class CRes : public QWidget {
	Q_OBJECT
public:
	CRes(map<string, string>&, QWidget *parent = 0);
	virtual ~CRes();
private:
	map<string, string> &paramsMap;
	//Res
	QLineEdit* delta_FF_set;
	QLineEdit* pion_axial_mass;
	QLineEdit* pion_C5A;
	//RES - DIS
	QLineEdit* spp_precision;
	QLineEdit* res_dis_cut;
private slots:
void delta_FF_setChanged(const QString&);
void pion_axial_massChanged(const QString&);
void pion_C5AChanged(const QString&);
void spp_precisionChanged(const QString&);
void res_dis_cutChanged(const QString&);
};

#endif /* CRES_H_ */
