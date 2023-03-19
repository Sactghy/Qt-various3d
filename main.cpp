#include "mainwindow.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    MainWindow w { &w };

    w.control->show();
    w.info->show();
    w.show();

    w.setFocus();


    return a.exec();
}
