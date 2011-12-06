/*
 * CCoh.h
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#ifndef CCOH_H_
#define CCOH_H_

#include <QtGui>
#include <map>
using namespace std;


class CCoh : public QWidget {
	Q_OBJECT
public:
	CCoh(map<string, string>&, QWidget *parent = 0);
	virtual ~CCoh();
private:
	map<string, string> &paramsMap;
	QLineEdit* coh_mass_correction;
private slots:
void coh_mass_correctionChanged(const QString&);
};

#endif /* CCOH_H_ */
