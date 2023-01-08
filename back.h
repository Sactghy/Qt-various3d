#ifndef BACK_H
#define BACK_H

#include <QDialog>

namespace Ui { class Back; }

class Back : public QDialog
{
    Q_OBJECT

public:
    explicit Back(QWidget *parent = nullptr);

    ~Back();

    Ui::Back *ui;

    bool changed{false};

private slots:


private:

};

#endif // BACK_H
