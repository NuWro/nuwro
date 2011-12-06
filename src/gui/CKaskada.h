/*
 * CKaskada.h
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#ifndef CKASKADA_H_
#define CKASKADA_H_

#include <QtGui>
#include <map>
using namespace std;

class CKaskada : public QLabel {
	Q_OBJECT
public:
	CKaskada(map<string, string>&, QWidget *parent = 0);
	virtual ~CKaskada();
private:
	map<string, string> &paramsMap;
	//nucleus_model       = 0
//	paramsMap["kaskada_debug"] = "";
	QCheckBox *kaskada_on;
	QCheckBox* pauli_blocking;
	QLineEdit* formation_zone;
	QComboBox* xsec;//model of cross sections in cascade: 0 - metropolis, 1 - oset
private slots:
void kaskada_onChanged();
void pauli_blockingChanged();
void formation_zoneChanged(const QString&);
void xsecChanged(const QString&);
};

#endif /* CKASKADA_H_ */
