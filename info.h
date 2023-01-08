#ifndef INFO_H
#define INFO_H

#include <QDialog>
#include "ui_info.h"

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

#endif // INFO_H
