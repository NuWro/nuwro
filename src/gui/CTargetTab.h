/*
 * CTargetTab.h
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#ifndef CTARGETTAB_H_
#define CTARGETTAB_H_

#include <QtGui>
#include <map>
#include "CElement.h"
using namespace std;


class CTargetTab : public QWidget{
	Q_OBJECT
private:
	map<string, string> &paramsMap;
	QComboBox *targetType;
    QComboBox *nucleus_model;
    vector<string> target_content;

    QLineEdit *geo_file;
    QLineEdit *geo_name;
    QLineEdit *geo_o;
    QLineEdit *geo_d;

	QGroupBox *targetTypeZeroGroup;
	QGroupBox *targetTypeOneGroup;
	QGroupBox *targetTypeTwoGroup;
public:
	CTargetTab(map<string, string> &, QWidget *parent = 0);
	virtual ~CTargetTab();
private slots:
	void targetTypeChanged(const QString&);
	void nucleus_modelChanged(const QString&);
	void geo_fileChanged(const QString& qs);
	void loadFileName();
	void geo_nameChanged(const QString& qs);
	void geo_oChanged(const QString& qs);
	void geo_dChanged(const QString& qs);

};


#endif /* CTARGETTAB_H_ */
