#ifndef CONTROL_H
#define CONTROL_H

#include <QDialog>
#include "ui_control.h"

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

    void on_pushRotate_clicked();

private:



};

#endif // CONTROL_H
