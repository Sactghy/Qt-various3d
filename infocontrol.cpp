#include "infocontrol.h"

InfoDial::InfoDial(QWidget *parent) : QDialog(parent), ui(new Ui::InfoDial)
{
    ui->setupUi(this);

    this->setGeometry( 1390, 575, this->height(), this->width() );

}

InfoDial::~InfoDial()
{
    delete ui;
}

#include "ui_control.h"
#include <math.h>


Dialog::Dialog(QWidget *parent) : QDialog(parent), ui(new Ui::Dialog)
{
    this->setGeometry( 1390, 150, this->height(), this->width() );

    ui->setupUi(this);//QLineEdit

}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::on_bQuit_clicked()
{
    QApplication::instance()->quit();
}

void Dialog::on_sld00_valueChanged(int value)
{
    sld00_v = ( 1.0 / 1000.0 ) * std::pow( ( static_cast<double>( value ) / 1000.0 ), ord  );
}

void Dialog::on_checkCur_clicked(bool checked)
{
     changed = true;
     cursor = checked;
}
