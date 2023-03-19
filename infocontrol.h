#ifndef INFOCONTROL_H
#define INFOCONTROL_H

#include <QDialog>
#include "ui_info.h"
#include "ui_control.h"

namespace Ui { class InfoDial; }

class InfoDial : public QDialog
{
    Q_OBJECT

public:

    explicit InfoDial(QWidget *parent = nullptr);

    ~InfoDial();

    Ui::InfoDial *ui;

    bool changed{false};


private slots:


private:


};


namespace Ui { class Dialog; }

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = nullptr);

    ~Dialog();

    Ui::Dialog *ui;

    bool changed{false}, cursor{false};

    double sld00_v{0.001}; double ord{10};

private slots:

    void on_bQuit_clicked();

    void on_sld00_valueChanged(int value);

    void on_checkCur_clicked(bool checked);

private:



};

#endif // INFOCONTROL_H
