#ifndef CMAINWINDOW_H_
#define CMAINWINDOW_H_
/*
 * CMainWindow.h
 *
 *  Created on: 31-05-2011
 *      Author: Marcin Boczulak
 */
#include<iostream>
#include<string.h>
#include<params.h>
#include<vector>
#include<map>
#include <nuwro.h>
#include <QtGui>

#include "CGeneralTab.h"
#include "CBeamTab.h"
#include "CTargetTab.h"
#include "CPhysicsTab.h"
using namespace std;



class CMainWindow : public QMainWindow
{
	Q_OBJECT
private:
	char **argv;

	map <string, string> paramsMap;
	QTabWidget *tab;
	string params;
	QPushButton *startButton;
	QPushButton *quitButton;
	QLayout *layout;
	QLayout *buttonLayout;
	QWidget *widget;
	QMenuBar *menu;


public:
	CMainWindow();
	virtual ~CMainWindow();
	void setDefaultParams();

public slots:
	void startNuwro();


};



#endif /* CMAINWINDOW_H_ */
