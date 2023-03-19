#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "infocontrol.h"
#include "geometry.h"

#include <sstream>
#include <algorithm>
#include <utility>
#include <thread>
#include <random>
#include <mutex>
#include <cassert>
#include <exception>
#include <fstream>

#include <QMainWindow>
#include <QString>
#include <QGraphicsScene>
#include <QGraphicsPixmapItem>
#include <QTimer>
#include <QKeyEvent>

class OD {
public: Object* o{}; double d{}; bool isB{};
        OD( Object* o0, double d0, bool i0 = false ) : o{o0}, d{d0}, isB{i0} {}
};

class ODD {
public: Object* o{}; double d{}; bool isB{}; Vec3d pnt{};
        ODD( Object* o0, double d0, bool i0 = false, Vec3d p0 = Vec3d{} ) : o{o0}, d{d0}, isB{i0}, pnt{p0} {}
};

class BSD {
public: bSphere* s{}; double d{}; bool isB{};
        BSD( bSphere* s0, double d0, bool i0 = false ) : s{s0}, d{d0}, isB{i0} {}
};

class DCO {
public: double distance{INFINITY}; Object* opnt{}; Vec3d colour{};
        DCO ( double d, Object* o, Vec3d c ) { distance = d; opnt = o; colour = c; }
};

class COO {
public: Object* opnt{}; Vec3d colour{}; double opq{};
        COO ( Vec3d c, Object* o, double o0 = 1.0 ) { opnt = o; colour = c; opq = o0; }
};

class CDO {
public: Vec3d colour{}; double opq{}, dist{};
        CDO ( Vec3d c, double d0 = INFINITY, double o0 = 1.0 ) { colour = c; dist = d0; opq = o0; }
};


class Rand
{
public:

unsigned int randi( unsigned int min, unsigned int max )
{	std::random_device rd; std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(min, max);
    return distrib(gen);  }

double randd( int min, int max )
{	std::random_device rd; std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(min, max);
    return distrib(gen);  }
};

class Timer
{
private:

    using clock_type = std::chrono::steady_clock;
    using second_type = std::chrono::duration< double, std::ratio<1,60> >;

    std::chrono::time_point<clock_type> m_beg;

public:
    Timer() : m_beg { clock_type::now() } { }

    void reset() { m_beg = clock_type::now(); }

    double elapsed() const
    { return std::chrono::duration_cast<second_type>( clock_type::now() - m_beg ).count(); }
};

class thGrd
{

std::thread t;

public: explicit thGrd( std::thread t_ ) : t( std::move(t_) )
        { if ( !t.joinable() ) throw std::logic_error("No thread!"); }

        ~thGrd() { t.join(); }

        thGrd( thGrd const& ) = delete;

        thGrd& operator=( thGrd const& ) = delete;
};

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    MainWindow(MainWindow *ptTh, QWidget *parent = nullptr);

    ~MainWindow();

    Dialog *control;
    InfoDial *info;

    double resc{0.001};

protected:

    void keyPressEvent(QKeyEvent *k)
    {
        kkey = k->key();
    };

    void keyReleaseEvent(QKeyEvent*)
    {
        kkey = 0;
    };

private slots:

    void startF();

    void endF();

    void rotateSph();

    void planeTxt( bool is );

    void blur();

    void createInbSpheres( Vec3d bpos );

private:

    Ui::MainWindow *ui;

    QWidget *pntThis;

    QGraphicsScene scene;

    uchar *dtX, *dtZ, *dtO, *dtJ; char *dtT;
    QImage *img0;
    QPixmap *pix0;
    QGraphicsPixmapItem *qGpixItem0;

    Timer timer;
    QTimer *qT0, *qT1;

    QCursor cursor;
    QPoint cPos{ cursor.pos() };

    const double pi180 { M_PI / 180 }, phi{ 1.618033988749 };

    int curX{ cPos.x() }, curY{ cPos.y() }, kcnt{}, kkey{};
    int rcnt{}, bcnt{}, gcnt{64};
    double looptime{}, degtorad;
    double qDegRadX{qDegreesToRadians(0.0)}, qDegRadY{qDegreesToRadians(0.0)}, qDegRadZ{qDegreesToRadians(180.0)};
    double sc0{1.0}, sc1{1.0}, sc2{1.0}, sc3{1.0};
    bool nomouse{false}, cursis{true};
    bool sb0{}, sb1{}, sb2{}, sb3{};
    bool ris{}, bis{}, gis{};


    double ksx[3][3] = { {-0.33, 0.00, 0.33}, {-0.69, 0.0, 0.69}, {-0.33, 0.00, 0.33} },
           ksy[3][3] = { {-0.33,-0.69,-0.33}, { 0.00, 0.0, 0.00}, { 0.33, 0.69, 0.33} };

    double deg{16}; bool yesno{false};

    Matrix44d cam, rmx;

    Object *pntd {nullptr};

    sSphere tObj;

    Rand rand;

    std::vector<bSphere*> rbSph;
    std::vector<oSphere*> randSph;
    std::vector<aTriangle*> tetra5;
    std::vector<Point*> pntarr;

    wSphere worldSphere { { }, 6999, { 0, 0, 0 }, 1 }, crs0 { { }, 1, { 0, 0, 99 }, 0.76 };

    bSphere bsp0 { {-500,100,-100}, 116, 0, nullptr }, bsp00 { bsp0.pos, 105, 0, &bsp0, oType::triangle };
    bSphere bsp1 { {-550,-50,-450}, 76, 0, nullptr, oType::sphere };
    bSphere bsp4 { {-600,-200,0}, 150, 0, nullptr, oType::plane };
    bSphere bsp5 { {-450,-150,-200}, 25, 1, nullptr, oType::line };
    bSphere bsp9 { {-250,-175,-300}, 120, 0, nullptr, oType::order2 };
    bSphere bspB { {-700,0,350}, 200, 0, nullptr, oType::order2 };
    bSphere bspT { {-700,-100,-250}, 100, 0, nullptr, oType::point };

    std::vector<bSphere*> bSpheres {  &bsp0, &bsp9 , &bspB, &bsp5, &bspT, &bsp4, &bsp1 }; //

    sSphere ss0 { bsp0.pos, 100, {66,77,36}, 0.35, &bsp0, 1 };

    Vec3d vsb{bsp1.pos.x+15, bsp1.pos.y+10, bsp1.pos.z+5},  vsq{bsp1.pos.x+10, bsp1.pos.y+25, bsp1.pos.z+35}, vsd{bsp1.pos.x+20, bsp1.pos.y+5, bsp1.pos.z-15},
          vsc{bsp1.pos.x-25, bsp1.pos.y+20, bsp1.pos.z+20}, vsw{bsp1.pos.x-15, bsp1.pos.y+15, bsp1.pos.z+20}, vsa{bsp1.pos.x-5, bsp1.pos.y+15, bsp1.pos.z+20},
          vsg{bsp1.pos.x+30, bsp1.pos.y-30, bsp1.pos.z+30}, vsr{bsp1.pos.x+20, bsp1.pos.y-10, bsp1.pos.z+25}, vsz{bsp1.pos.x+10, bsp1.pos.y-15, bsp1.pos.z+15},
          vse{bsp1.pos.x-20, bsp1.pos.y-10, bsp1.pos.z-20}, vsy{bsp1.pos.x-10, bsp1.pos.y-30, bsp1.pos.z-20}, vsx{bsp1.pos.x-15, bsp1.pos.y-25, bsp1.pos.z-25},
          vsi{bsp1.pos.x-35, bsp1.pos.y+10, bsp1.pos.z-25}, vsf{bsp1.pos.x-15, bsp1.pos.y+30, bsp1.pos.z-35}, vsm{bsp1.pos.x-15, bsp1.pos.y+25, bsp1.pos.z-25},
          vsn{bsp1.pos.x+35, bsp1.pos.y-20, bsp1.pos.z-30}, vss{bsp1.pos.x+35, bsp1.pos.y-20, bsp1.pos.z-10}, vsv{bsp1.pos.x+25, bsp1.pos.y-20, bsp1.pos.z-20};

    oSphere ss1b  { vsb, 23, {155,69,43}, 1, &bsp1 },   ss1q  { vsq, 13, {135,76,45}, 1, &bsp1 },   ss1d  { vsd, 16, {115,69,55}, 1, &bsp1 },
            ss1c  { vsc, 19, {77,55,99}, 1, &bsp1 },    ss1w  { vsw, 11, {69,54,99}, 1, &bsp1 },    ss1a  { vsa, 19, {44,54,111}, 1, &bsp1 },
            ss1g  { vsg, 16, {99,84,45}, 1, &bsp1 },    ss1r  { vsr, 13, {96,88,55}, 1, &bsp1 },    ss1z  { vsz, 13, {99,66,55}, 1, &bsp1 },
            ss1e  { vse, 19, {44,64,95}, 1, &bsp1 },    ss1y  { vsy, 11, {39,63,93}, 1, &bsp1 },    ss1x  { vsx, 11, {42,31,111}, 1, &bsp1 },
            ss1i  { vsi, 19, {64,77,99}, 1, &bsp1 },    ss1f  { vsf, 13, {66,77,99}, 1, &bsp1 },    ss1m  { vsm, 13, {69,42,99}, 1, &bsp1 },
            ss1n  { vsi, 16, {78,55,111}, 1, &bsp1 },   ss1s  { vss, 11, {69,54,111}, 1, &bsp1 },   ss1v  { vsv, 16, {63,76,111}, 1, &bsp1 };

    Vec3d t0 { 1.0, 1.0, 1.0 },
          t1 { 1.0,-1.0,-1.0 },
          t2 {-1.0, 1.0,-1.0 },
          t3 {-1.0,-1.0, 1.0 };

    Plane pln0 { bsp4.pos, { 247, 169, 0 }, { 34, 166, 99 }, 1, &bsp4 };

    Line  pOx  { bsp5.pos, { bsp5.pos.x+25, bsp5.pos.y, bsp5.pos.z }, { 111,136,99 }, 1, &bsp5 },
          pOy  { bsp5.pos, { bsp5.pos.x, bsp5.pos.y+25, bsp5.pos.z }, { 134,166,99 }, 1, &bsp5 },
          pOz  { bsp5.pos, { bsp5.pos.x, bsp5.pos.y, bsp5.pos.z+25 }, { 144,196,99 }, 1, &bsp5 },
          mOx  { bsp5.pos, { bsp5.pos.x-25, bsp5.pos.y, bsp5.pos.z }, { 39,69,136 }, 1, &bsp5 },
          mOy  { bsp5.pos, { bsp5.pos.x, bsp5.pos.y-25, bsp5.pos.z }, { 69,43,169 }, 1, &bsp5 },
          mOz  { bsp5.pos, { bsp5.pos.x, bsp5.pos.y, bsp5.pos.z-25 }, { 39,16,196 }, 1, &bsp5 };

    Ellipse l20 { bsp9.pos, &bsp9, 0.63, { 99.0, 159.0, 135.0 }, 500, 3 };
    Ellipse l21 { bsp9.pos, &bsp9, 0.48, { 48.0, 135.0, 144.0 }, 800, 1 };
    Ellipse l22 { bsp9.pos, &bsp9, 0.69, { 99.0, 105.0, 165.0 }, 1000, 4 };
    Ellipse l23 { bsp9.pos, &bsp9, 0.48, { 169.0, 48.0, 135.0 }, 1400, 2 };

    Hyper1 h1x { bspB.pos, &bspB, { 65.0, 95.0, 165.0 } };
    Hyper1 h2x { bspB.pos, &bspB, { 75.0, 65.0, 195.0 } };
    Hyper2 h2y { bspB.pos, &bspB, { 95.0, 135.0, 95.0 } };

    Vec3d from{ -900.0, 75.0, -500.0 }, to{ bspT.pos };

};
#endif // MAINWINDOW_H




