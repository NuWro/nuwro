/*
 * CGeneralTab.h
 *
 *  Created on: 19-06-2011
 *      Author: Marcin Boczulak
 */

#ifndef CGENERALTAB_H_
#define CGENERALTAB_H_

#include <QtGui>
#include <map>
using namespace std;

class CGeneralTab: public QWidget {
Q_OBJECT
private:
	map<string, string> &paramsMap;
	QLineEdit *numberOfEvents;
	QLineEdit *numberOfTestEvents;
	QLineEdit *randomSeed;
	QCheckBox *dyn_qel_cc; // Quasi elastic charged current
	QCheckBox *dyn_qel_nc; // Quasi elastic neutral current
	QCheckBox *dyn_res_cc; // Resonant charged current
	QCheckBox *dyn_res_nc; // Resonant neutral current
	QCheckBox *dyn_dis_cc; // Deep inelastic charged current
	QCheckBox *dyn_dis_nc; // Deep inelastic neutral current
	QCheckBox *dyn_coh_cc; // Coherent charged current
	QCheckBox *dyn_coh_nc; // Coherent neutral current



public:
	CGeneralTab(map<string, string>&, QWidget *parent = 0);
	virtual ~CGeneralTab();

private slots:
	void numberOfEventsChanged(const QString&);
	void numberOfTestEventsChanged(const QString&);
	void randomSeedChanged(const QString&);
	void setDynQelCc();
	void setDynQelNc();
	void setDynResCc();
	void setDynResNc();
	void setDynDisCc();
	void setDynDisNc();
	void setDynCohCc();
	void setDynCohNc();

};

#endif /* CGENERALTAB_H_ */
