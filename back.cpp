#include "back.h"
#include "ui_back.h"


Back::Back(QWidget *parent) : QDialog(parent), ui(new Ui::Back)
{
    ui->setupUi(this);

    this->setGeometry( 1, 0, this->height(), this->width() );


}

Back::~Back()
{
    delete ui;
}

