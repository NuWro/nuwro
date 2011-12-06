/*
 * CGel.h
 *
 *  Created on: 29-06-2011
 *      Author: boczus
 */

#ifndef CQEL_H_
#define CQEL_H_

#include <QtGui>
#include <map>
using namespace std;

class CQel : public QWidget{
	Q_OBJECT
public:
	CQel(map<string, string>&, QWidget *parent = 0);
	virtual ~CQel();
private:
	map<string, string> &paramsMap;
	QSpinBox *qel_new;  // 0 - old qel implementaions, 1 - new qel implementation
	QLineEdit *qel_cc_axial_mass; //[MeV] axial mass
	QLineEdit *qel_nc_axial_mass; //[MeV] axial mass
	QSpinBox *qel_cc_vector_ff_set;//electromagnetic FF: 1 -> dipole; 2 -> BBBA;
	QSpinBox *qel_cc_axial_ff_set; // 1=dipole, 2,3,4... =  2-fold parabolic modification of axial FF
	QSpinBox *qel_strange;//0-1
	QLineEdit *delta_s;
	QLineEdit *qel_s_axial_mass;
	QSpinBox *flux_correction;//0-1
	QSpinBox *sf_method;//0-2
	QSpinBox *sf_form_factors;// 0-3
	QSpinBox *cc_smoothing;//0-1
//	QSpinBox *qel_kinematics;// 0-4


private slots:
void qel_newChanged(const QString&);
void qel_cc_axial_massChanged(const QString&);
void qel_nc_axial_massChanged(const QString&);
void qel_cc_vector_ff_setChanged(const QString&);
void qel_cc_axial_ff_setChanged(const QString&);
void qel_strangeChanged(const QString&);
void delta_sChanged(const QString&);
void qel_s_axial_massChanged(const QString&);
void flux_correctionChanged(const QString&);
void sf_methodChanged(const QString&);
void sf_form_factorsChanged(const QString&);
void cc_smoothingChanged(const QString&);


};

#endif /* CQEL_H_ */
