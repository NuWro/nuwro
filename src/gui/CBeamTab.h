/*
 * CBeamTab.h
 *
 *  Created on: 21-06-2011
 *      Author: boczus
 */

#ifndef CBEAMTAB_H_
#define CBEAMTAB_H_

#include <QtGui>
#include <map>

using namespace std;


class CBeamTab : public QWidget{
	Q_OBJECT
private:
	map<string, string> &paramsMap;
	QComboBox *beamType;
	QLineEdit *beamDirection;
	QLineEdit *beamParticle;
	QLineEdit *beamEnergy;
	QLineEdit *beamContent;
	QLineEdit *beamFileName;

	QVBoxLayout *vBoxBeamTypeZero;
	QVBoxLayout *vBoxBeamTypeOne;
	QVBoxLayout *vBoxBeamTypeTwo;
	QVBoxLayout *vBoxBeamTypeThree;
	QSpacerItem *vSpacer;

	QGroupBox *beamTypeZeroGroup;
	QGroupBox *beamTypeOneGroup;
	QGroupBox *beamTypeTwoGroup;
	QGroupBox *beamTypeThreeGroup;

public:
	CBeamTab(map<string, string> &, QWidget *parent = 0);
	virtual ~CBeamTab();

private slots:
	void beamTypeChanged(const QString&);
	void beamDirectionChanged(const QString&);
	void beamParticleChanged(const QString&);
	void beamEnergyChanged(const QString&);
	void beamContentChanged(const QString&);
	void beamFileNameChanged(const QString&);
	void loadFileName();
};

#endif /* CBEAMTAB_H_ */
