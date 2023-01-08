#include "info.h"


InfoDial::InfoDial(QWidget *parent) : QDialog(parent), ui(new Ui::InfoDial)
{
    ui->setupUi(this);

    this->setGeometry( 1390, 575, this->height(), this->width() );

}

InfoDial::~InfoDial()
{
    delete ui;
}

