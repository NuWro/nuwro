
#include "nuwro.h"
#include "CMainWindow.h"
#include <QtGui/QApplication>
#include <QTableView>
#include <iostream>
#include <QtGui>
//#include <QCoreApplication>
//#include <QtNuwro.h>

int main(int argc,char *argv[])
{
    QApplication app(argc, argv);
    CMainWindow window;
    window.show();
    return app.exec();
}
