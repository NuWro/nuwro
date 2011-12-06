/*
 * CElement.h
 *
 *  Created on: 26-06-2011
 *      Author: boczus
 */

#ifndef CELEMENT_H_
#define CELEMENT_H_

#include <QtGui>
#include <map>
#include <qabstractspinbox.h>
using namespace std;

class CElement : public QWidget
{
	Q_OBJECT
public:
    CElement(map<string,string>& , QWidget *parent = 0);
    CElement(map<string,string>& , vector<string>, QWidget *parent = 0);
    virtual ~CElement();
private:
	map<string, string> &paramsMap;
	QSpinBox *nucleus_p;
	QSpinBox *nucleus_n;
	QLineEdit *nucleus_E_b;
	QLineEdit *nucleus_kf;
	QComboBox *nucleus_target;

	QSpinBox *nucleus_p2;
	QSpinBox *nucleus_n2;
	QSpinBox *count;
	QSpinBox *count2;
	QLineEdit *nucleus_E_b2;
	QLineEdit *nucleus_kf2;
	QComboBox *nucleus_target2;

	vector<string> target_content;

	QVBoxLayout *vBoxTargetTypeZero;

    void changeParamTargetContent();

private slots:
	void nucleus_pChanged(const QString&);
	void nucleus_nChanged(const QString&);
	void nucleus_E_bChanged(const QString&);
	void nucleus_kfChanged(const QString&);
	void nucleus_targetChanged(const QString&);

	void nucleus_pTargetContentChanged(const QString&);
	void nucleus_nTargetContentChanged(const QString&);
	void countTargetContentChanged(const QString&);
	void nucleus_E_bTargetContentChanged(const QString&);
	void nucleus_kfTargetContentChanged(const QString&);
	void nucleus_targetTargetContentChanged(const QString&);
	void nucleus_p2TargetContentChanged(const QString&);
	void nucleus_n2TargetContentChanged(const QString&);
	void count2TargetContentChanged(const QString&);
	void nucleus_E_b2TargetContentChanged(const QString&);
	void nucleus_kf2TargetContentChanged(const QString&);
	void nucleus_target2TargetContentChanged(const QString&);
};

#endif /* CELEMENT_H_ */
