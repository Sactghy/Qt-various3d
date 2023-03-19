#include "mainwindow.h"
#include "ui_mainwindow.h"

Vec3d fvec, orig;
bool is_pln0;
double tst{1};
Object *transObj;

MainWindow::MainWindow( MainWindow *ptTh, QWidget *parent ) : QMainWindow(parent), ui(new Ui::MainWindow), pntThis(ptTh)
{
    num_threads = std::thread::hardware_concurrency();

    num_threads = ( num_threads >= 12 ) ? 12 : num_threads;

    bsp5.opq = 0.05;

    tObj.opq = 1.0; transObj = dynamic_cast<class Object*>(&tObj);

    // 5 TETRA

    Vec3d t00{},t01{},t02{},t03{},te0{t0},te1{t1},te2{t2},te3{t3};

    t0 = te0 * 55.0; t1 = te1 * 55.0; t2 = te2 * 55.0; t3 = te3 * 55.0;

    uint icnt = {};

    do { Vec3d c1{55,135,99}, c2{88,99,111}, c3{99,66,177}, c4{155,55,99};

    if ( icnt == 1 || icnt == 3 ) { c2 = {155,55,99}; c4 = {99,66,177}; c3 = {55,135,99}; c1 = {88,99,111}; }
    else if ( icnt == 2 ) { c2 = {55,135,99}; c4 = {88,99,111}; c3 = {155,55,99}; c1 = {99,66,177}; }

    t00 = t0 + bsp00.pos; t01 = t1 + bsp00.pos; t02 = t2 + bsp00.pos; t03 = t3 + bsp00.pos;

    aTriangle* t5t0 = new aTriangle{ t00, t01, t02, 1, c1, &bsp00 }; tetra5.push_back(t5t0);
    aTriangle* t5t1 = new aTriangle{ t00, t02, t03, 1, c2, &bsp00 }; tetra5.push_back(t5t1);
    aTriangle* t5t2 = new aTriangle{ t01, t00, t03, 1, c3, &bsp00 }; tetra5.push_back(t5t2);
    aTriangle* t5t3 = new aTriangle{ t01, t02, t03, 1, c4, &bsp00 }; tetra5.push_back(t5t3);

    icnt++; } while ( icnt < 5 );

    rmx.m[0][0] = ( phi - 1.0 ) / 2.0;    rmx.m[0][1] = phi / 2.0;              rmx.m[0][2] =  1.0 / 2.0;
    rmx.m[1][0] = ( -1.0 * phi ) / 2.0;   rmx.m[1][1] = 1.0 / 2.0;              rmx.m[1][2] = ( 1.0 - phi ) /2.0;
    rmx.m[2][0] = -1.0 / 2.0;             rmx.m[2][1] = ( 1.0 - phi ) / 2.0;    rmx.m[2][2] = phi /2.0;

    icnt = {};
    for ( auto t5 : tetra5 ) { if ( icnt > 5 && icnt % 4 == 0 ) rmx = rmx * rmx; if ( icnt >= 4 )
        { t0 = t5->a - bsp00.pos; t1 = t5->b - bsp00.pos; t2 = t5->c - bsp00.pos;
          rmx.multVecMatrix(t0,te0); rmx.multVecMatrix(t1,te1); rmx.multVecMatrix(t2,te2);
          t5->a = te0 + bsp00.pos; t5->b = te1 + bsp00.pos; t5->c = te2 + bsp00.pos; }
    icnt++; }

    // RECT

    std::srand( static_cast<int>( rand.randi( 1, 16634 ) ) ); std::rand();

    dtO = new uchar[494*339*3]{};
    dtZ = new uchar[494*339*3]{};
    dtJ = new uchar[494*339*3]{}; pln0.dtZ = dtJ;

   for ( int x = 0; x <= 494; x++ ) { for ( int y = 0; y <= 338; y++ ) {

        int p = ( ( y * 494 ) + x ) * 3;

     dtO[p+0] = std::rand();  dtO[p+1] = std::rand();  dtO[p+2] = std::rand();

        dtZ[p+0] = dtO[p+0];  dtZ[p+1] = dtO[p+1];  dtZ[p+2] = dtO[p+2]; } }

    for ( int n = 4; n != 0; n-- ) planeTxt( true );

    for ( int x = 0; x <= 494; x++ ) { for ( int y = 0; y <= 338; y++ ) {

         int p = ( ( y * 494 ) + x ) * 3;

      dtO[p+0] += 96; dtO[p+1] += 128; dtO[p+2] += 64; } }

    // DOTS

    std::vector<Vec3d> dircol;

    Matrix44d a;
    double xx{16},yy{16},zz{16};

    for ( uint nn = 0; nn < 8; nn++ ) { Vec3d col{66,199,133}, pos{1,1,1}, ppp{}; double size{1};
                                        double deg{ rand.randd( 1400, 1574 ) }; bool yynn{false};
                                        while( static_cast<int>(deg) % 25 != 0 ) deg++;
                                        if ( deg > 1500 ) yynn = true; deg /= 100;

    for ( uint ii = 0; ii < 39; ii++ ) {

        double qDegRadX = qDegreesToRadians(xx);
        double qDegRadY = qDegreesToRadians(yy);
        double qDegRadZ = qDegreesToRadians(zz);

        a.m[0][0] = qCos(qDegRadX) * qCos(-qDegRadY); a.m[0][1] = ( qCos(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) - ( qSin(qDegRadX) * qCos(qDegRadZ) ); a.m[0][2] = ( qCos(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) + ( qSin(qDegRadX) * qSin(qDegRadZ) );
        a.m[1][0] = qSin(qDegRadX) * qCos(qDegRadY);  a.m[1][1] = ( qSin(-qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) + ( qCos(-qDegRadX) * qCos(qDegRadZ) ); a.m[1][2] = ( qSin(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) - ( qCos(qDegRadX) * qSin(qDegRadZ) );
        a.m[2][0] = -1.0 * qSin(qDegRadY);            a.m[2][1] = qCos(qDegRadY) * qSin(qDegRadZ);                                                            a.m[2][2] = qCos(qDegRadY) * qCos(qDegRadZ);

        Point* p0 = new Point { ppp + bspT.pos, { col.x + rand.randd( -19, 31 ),
                                                  col.y + rand.randd( -31, 54 ),
                                                  col.z + rand.randd( -31, 33 ) },
                                rand.randd( 69, 100 ) / 100.0, nullptr, size };

        if ( ii > 19 ) size *= 0.87;
        p0->deg = deg; p0->yn = yynn;

        pntarr.push_back(p0);

        pos = pos * 1.087;
        a.multDirMatrix( pos, ppp );
        pos = ppp;
        switch ( nn ) {
            case 0 : break;
            case 1 : ppp.x *= -1; break;
            case 2 : ppp.y *= -1; break;
            case 3 : ppp.z *= -1; break;
            case 4 : ppp.x *= -1; ppp.z *= -1; break;
            case 5 : ppp.y *= -1; ppp.z *= -1; break;
            case 6 : ppp.x *= -1; ppp.y *= -1; break;
            case 7 : ppp.x *= -1; ppp.y *= -1; ppp.z *= -1; break;
        }

    } }

    Vec3d bpos { bspT.pos }; icnt = {16};

    while ( icnt != 2 ) {

        createInbSpheres( bpos );

    switch ( icnt ) {

    case 16 : bpos = bspT.pos + 13; break;
    case 15 : bpos = bspT.pos - 13; break;
    case 14 : bpos = bspT.pos + 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; break;
    case 13 : bpos = bspT.pos - 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; break;

    case 12 : bpos = bspT.pos - 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 11 : bpos = bspT.pos + 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 10 : bpos = bspT.pos - 13; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 9  : bpos = bspT.pos + 13; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;

    case 8 : bpos = bspT.pos; bpos.z += 23; break;
    case 7 : bpos = bspT.pos; bpos.z -= 23; break;
    case 6 : bpos = bspT.pos; bpos.y += 23; break;
    case 5 : bpos = bspT.pos; bpos.y -= 23; break;
    case 4 : bpos = bspT.pos; bpos.x += 23; break;
    case 3 : bpos = bspT.pos; bpos.x -= 23; break;
    }

    icnt--; } createInbSpheres( bpos );

    for ( auto p0 : pntarr ) { bool isin { false };

        for ( auto b0 : bspT.bSpheres ) {

            if ( ( p0->pos - b0->pos ) . length() < b0->rad - 1 ) {

                for ( auto b1 : b0->bSpheres ) {

                    if ( ( p0->pos - b1->pos ) . length() < b1->rad - 3.1 ) {

                       if ( !p0->owner ) { b1->owns.push_back(p0); isin = true; p0->owner = b1; }

                    } if ( isin ) break;

                } if ( !isin && !p0->owner ) { std::cout << "missed <- \n"; b0->owns.push_back(p0);
                       p0->owner = b0; p0->col = { 255, 255, 255 }; isin = true; }

            }

        }

    }

    for ( auto b0 : bspT.bSpheres ) {

        for ( std::vector<bSphere*>::iterator it = b0->bSpheres.begin(); it != b0->bSpheres.end();) {

            if ( it.operator*()->owns.empty() ) { b0->bSpheres.erase( it ); }
            else {  Vec3d newpos{}; double ocnt{};
                    for ( auto o : it.operator*()->owns ) {
                    newpos = newpos + ( o->pos - it.operator*()->pos );
                    ocnt++;
                    } newpos = { newpos.x / ocnt, newpos.y / ocnt, newpos.z / ocnt };
                    it.operator*()->pos = newpos + it.operator*()->pos;
                    double res{};
                    for ( auto o : it.operator*()->owns ) {
                    if ( res < ( o->pos - it.operator*()->pos ) . length() + 1  )
                       { res = ( o->pos - it.operator*()->pos ) . length() + 1; }
                    } it.operator*()->rad = res; it.operator*()->rad2 = res * res;
                    if ( res > 10.5 ) std::cout << res << "\n";
                    it++; }

        } if ( !b0->owns.empty() && b0->owns.size() > 1 ) { std::cout << "not empty - \n"; }

    }  for ( std::vector<bSphere*>::iterator it = bspT.bSpheres.begin(); it != bspT.bSpheres.end();)
           { if ( it.operator*()->bSpheres.empty() ) bspT.bSpheres.erase(it);
             else { Vec3d newpos{  }; double ocnt{};
                    for ( auto b : it.operator*()->bSpheres ) {
                    newpos = newpos + ( b->pos - it.operator*()->pos ); ocnt++;
                    } newpos = { newpos.x / ocnt, newpos.y / ocnt, newpos.z / ocnt };
                    it.operator*()->pos = newpos + it.operator*()->pos;
                    double res{};
                    for ( auto b : it.operator*()->bSpheres ) {
                        if ( res < ( b->pos - it.operator*()->pos ) . length() + b->rad + 1 )
                           { res = ( b->pos - it.operator*()->pos ) . length() + b->rad + 1; }
                    } it.operator*()->rad = res; it.operator*()->rad2 = res * res;
                    it++; } }


    // O2

    h1x.s2i =3.9; h1x.s2j = -2.4; h1x.s2k = -3.6;
    h1x.s1i = 2.4; h1x.s1j = 3.3; h1x.s1k = 1.9;
    h1x.sci = 3.6; h1x.scj = 2.4; h1x.s1k = -2.4;
    h1x.sepr = true; h1x.opq = 0.69;

    h2x.s2i = 3.9; h2x.s2j = -2.4; h2x.s2k = -3.6;
    h2x.s1i = 2.4; h2x.s1j = 3.3; h2x.s1k = 1.9;
    h2x.sci = 9.9; h2x.scj = 2.4; h2x.s1k = -2.4;
    h2x.sepr = true; h2x.opq = 0.36;

    h2y.sci = -1.0; h2y.sck = -1.0; h2y.s1j = 4.0; h2y.s2i = 2.0;
    h2y.sect = true;

    // INIT

    ui->setupUi(this);

    control = new Dialog{pntThis}; info = new InfoDial{pntThis};

    scene.setSceneRect( -390, -390, 779, 779 );
    ui->graphicsView->setScene( &scene );
    ui->graphicsView->setFocusProxy( pntThis );

    dtX = new (std::nothrow) uchar[780*780*4] { 0 };
    for ( int i = 0; i <= 780*780*4; i += 4 ) { dtX[i+3]=255; }
    img0 = new QImage( dtX, 780, 780, QImage::Format_ARGB32 );
    pix0 = new QPixmap;
    pix0->convertFromImage( *img0, 0 );
    qGpixItem0 = scene.addPixmap( *pix0 );
    qGpixItem0->setPos( -390, -390 );

    dtT = new (std::nothrow) char[1360*1360*3]{};
    std::ifstream imgF("/zxc/qt/my/O/txh2.ppm", std::ios::in);
    if (!imgF.is_open()) std::cout << "open file: txh2.ppm error...\n";
    imgF.seekg( 0x40 ); while (imgF) { imgF.read( &dtT[0], (1360*1360*3) ); } imgF.close();

    qT0 = new QTimer(this); qT1 = new QTimer(this);
    connect(qT0, SIGNAL(timeout()), this, SLOT(startF()));
    connect(qT1, SIGNAL(timeout()), this, SLOT(endF()));
    qT0->start(63); qT1->start(0);

}

MainWindow::~MainWindow()
{
    delete ui; delete control; delete info;
    delete[] dtX; delete img0; delete pix0; delete qGpixItem0;
    delete[] dtZ; delete[] dtO; delete[] dtJ; delete[] dtT;

    for ( auto rs : randSph ) delete rs;

    for ( auto rb : rbSph ) delete rb;

    for ( auto t5 : tetra5 ) delete t5;

    for ( auto p0 : pntarr ) delete p0;

}


void MainWindow::startF()
{
    timer.reset();

    if ( kcnt == 0 ) switch ( kkey ) {

    case 47 : to = bsp00.pos; kcnt = 3; break;
    case 32 : nomouse = !nomouse; kcnt = 3; break;

    } else { kcnt--; }

    cPos = cursor.pos();

    int dx { curX - cPos.x() }, dy { ( curY - cPos.y() ) }; curX = cPos.x(); curY = cPos.y();

    dy = ( dy < 5 && dy > -5 ) ? 0 : static_cast<double>(dy) / 5.0;
    dx = ( dx < 5 && dx > -5 ) ? 0 : static_cast<double>(dx) / 5.0;

    Vec3d atpos = ( to - from ) . normalize();

    if ( dy != 0 && !nomouse ) { atpos.y += dy * pi180; to = from + atpos; }

    if ( dx != 0 && !nomouse ) { double xsin { sin( dx * pi180 ) }, xcos { cos( dx * pi180 ) },
                                 tx { ( atpos.x * xcos ) + ( atpos.z * xsin ) },
                                 tz { ( atpos.z * xcos ) - ( atpos.x * xsin ) };
                                 atpos.x = tx; atpos.z = tz; to = from + atpos; }

    Vec3d forward = ( from - to ) . normalize(),
          right = ( Vec3d( 0, 1, 0 ) . normalize() ) . cross( forward ),
          up = forward . cross( right );

    switch ( kkey )
    {
       case 16777235 : from = from - forward * 9; to = to - forward * 9; break;
       case 16777237 : from = from + forward * 9; to = to + forward * 9; break;

       case 16777236 : from = from + right * 9; to = to + right * 9; break;
       case 16777234 : from = from - right * 9; to = to - right * 9; break;
    }

    cam[0][0] = right.x;   cam[0][1] = right.y;   cam[0][2] = right.z;
    cam[1][0] = up.x;      cam[1][1] = up.y;      cam[1][2] = up.z;
    cam[2][0] = forward.x; cam[2][1] = forward.y; cam[2][2] = forward.z;
    cam[3][0] = from.x;    cam[3][1] = from.y;    cam[3][2] = from.z;

    Vec3d dest; std::vector<BSD> currentBVols {}; is_pln0 = false;

    fvec = to - from; orig = from;

    for ( auto b : bSpheres ) {  Vec3d tmp{ b->pos - from }; b->vis = false;

          double dp = fvec.dot( tmp ), curlen{ tmp.length() };
          bool isinside { curlen < b->rad };

          if ( ( dp > 0 ) || ( isinside && dp < 0 ) )
             { currentBVols.push_back( BSD { b, curlen - b->rad, isinside } ); b->vis = true; }
    }

    struct { bool operator()( BSD a, BSD b ) const { return a.d < b.d; } } sortBSD;
    std::sort( currentBVols.begin(), currentBVols.end(), sortBSD );

    auto processLine = [ & ] ( int y ) { double pX, pY; Vec3d dir;

      for ( int x = 0; x <=780; x++ ) { int p = ( ( 780 * y ) + x ) * 4; bool cursor {false};

       pX = ( 2 * ( ( x + 0.5 ) / 780 ) - 1 ); pY = ( 1 - 2 * ( ( y + 0.5 ) / 780 ) );

       cam.multDirMatrix( Vec3d( pX, pY, -1 ), dir ); dir.normalize();

       std::vector<COO> ray_res;

       if ( cursis && x >= 389 && x <= 391 && y >= 389 && y <= 391 )
       { cursor = true; ray_res.push_back( COO { crs0.col, &crs0, 1 - ( abs( 390.0 - x ) * 0.2 ) - ( abs( 390.0 - y ) * 0.2 ) } ); }

       double is1opq{};

         struct it_BVols { double dist1{INFINITY}, dist2{INFINITY}, prevLen{};

            void operator()( std::vector<COO> &ray_res,
                             std::vector<BSD> *currentBVols,
                             const Vec3d &dir,
                             double &is1opq )

            {   for ( auto b : *currentBVols ) {

                 double newLen {b.d}; Vec3d newOrig {orig};

                   do {

                      Vec3d futureOrig = newOrig + dir * newLen;
                      prevLen = newLen;
                      newOrig = futureOrig;
                      newLen = ( futureOrig - b.s->pos ).length() - b.s->rad;

                      if ( newLen < 0.1 )

                      { if ( b.s->isvis && !b.isB ) ray_res.push_back( COO { b.s->col, b.s } );

                        struct it_Objects { double dist1{INFINITY}, dist2{INFINITY}, prevLen{}; Vec3d res_col;

                         std::vector<OD> currentOwns {};
                         std::vector<ODD> currentOwnsO {};
                         std::vector<DCO> pre_rrs {};
                         std::vector<CDO> pre_cdo {};

                         void operator()( std::vector<COO> &ray_res,
                                          bSphere *bSphere,
                                          const Vec3d &dir,
                                          double &is1opq )
                        {

                  switch ( bSphere->otype ) {

                  case oType::none : {

                           for ( auto o : bSphere->owns ) { Vec3d tmp{ o->pos - orig };

                              double dp = fvec.dot( tmp ), curlen{ tmp.length() }, dist{ curlen - o->rad };

                              if ( ( dp > 0 ) || ( dist < 0 && dp < 0 ) ) currentOwns.push_back( OD { o, dist, false } );  }

                            for ( auto b : bSphere->bSpheres ) { Vec3d tmp{ b->pos - orig };

                              double dp = fvec.dot( tmp ), curlen{ tmp.length() }, dist{ curlen - b->rad };

                              if ( ( dp > 0 ) || ( dist < 0 && dp < 0  ) ) currentOwns.push_back( OD { b, dist, true } );  }


                            struct { bool operator()( OD a, OD b ) const { return a.d < b.d; } } sortOD;
                            std::sort( currentOwns.begin(), currentOwns.end(), sortOD );


                            for ( auto &o : currentOwns ) {

                                Vec3d newOrig{orig}; double newLen {o.d}, icnt{0};

                                do {

                                   Vec3d futureOrigCur = newOrig + dir * newLen;
                                   prevLen = newLen;
                                   newOrig = futureOrigCur;
                                   newLen = ( futureOrigCur - o.o->pos ).length() - o.o->rad;

                                  if ( newLen < 0.1  )
                                   {
                                       if ( o.isB ) {

                                           if ( o.o->isvis ) ray_res.push_back( COO { o.o->col, o.o } );
                                           it_Objects it; it( ray_res, dynamic_cast<class bSphere*>(o.o), dir, is1opq );
                                           if ( o.o->isvis ) ray_res.push_back( COO { o.o->col, o.o } );

                                       }

                                       else {

                                           Vec3d pnt = ( futureOrigCur - o.o->pos ) . normalize();
                                           res_col = { o.o->col.x + ( pnt.x / M_PI ) * 100,
                                                       o.o->col.y - ( pnt.y / M_PI ) * 100,
                                                       o.o->col.z + ( pnt.z / M_PI ) * 100 };

                                            ray_res.push_back( COO { res_col, o.o } );
                                            ray_res.push_back( COO { res_col, o.o } ); is1opq += o.o->opq;

                                       }

                                       newLen = INFINITY;

                                   } else if ( o.o->blur ) { icnt++;
                                          if ( icnt > 5 && icnt < 25 ) ray_res.push_back( COO { o.o->col, o.o, icnt/(icnt*10) } ); }

                                   } while ( newLen < prevLen ); if ( is1opq > 0.96 ) break;

                             } break; }

                  case oType::triangle : { int icnt{0};

                                for ( auto o : bSphere->owns) {

                                    if ( o->intersect( orig, dir, dist1, dist2, res_col ) )
                                    { pre_rrs.push_back( DCO { dist1, o, res_col } ); icnt++; /* if ( o->opq == 1 ) */ is1opq += o->opq; }

                                }

                                if ( icnt > 0 ) { if ( icnt > 1 ) {
                                struct { bool operator()( DCO a, DCO b ) const { return a.distance < b.distance; } } sortDCO;
                                std::sort( pre_rrs.begin(), pre_rrs.end(), sortDCO ); }
                                ray_res.push_back( COO { pre_rrs[0].colour, pre_rrs[0].opnt, pre_rrs[0].opnt->opq } );
                                }

                                break; }

                  case oType::sphere : {

                              for ( auto *o : bSphere->owns ) {

                                 double  curlen{ ( o->pos - orig ).length() }, dist{ curlen - o->rad };

                                 currentOwnsO.push_back( ODD { o, dist, false } ); }

                              for ( auto &o : currentOwnsO ) {

                                  Vec3d newOrig{orig}; double newLen { o.d };

                                  do {

                                      Vec3d futureOrigCur = newOrig + dir * newLen;
                                      prevLen = newLen;
                                      newOrig = futureOrigCur;
                                      newLen = ( futureOrigCur - o.o->pos ) . length() - o.o->rad;

                                      if ( newLen < 0.01 ) {

                                          o.pnt = futureOrigCur;
                                          o.d = ( futureOrigCur - orig ) . length();
                                          newLen = INFINITY;
                                          o.isB = true;
                                          is1opq = 1.0;
                                      }

                                  } while ( newLen < prevLen );

                                }

                                struct { bool operator()( ODD a, ODD b ) const { return a.d < b.d; } } sortODD;
                                std::sort( currentOwnsO.begin(), currentOwnsO.end(), sortODD );

                              for ( auto &o : currentOwnsO ) {

                                  if ( o.isB  ) {

                                      Vec3d pnt = ( o.pnt - o.o->pos ) . normalize();
                                      res_col = { o.o->col.x + atan ( pnt.x / M_PI ) * 100,
                                                  o.o->col.y - asinh ( pnt.y / M_PI ) * 100,
                                                  o.o->col.z + cosh ( pnt.z  / M_PI ) * 100 };

                                    ray_res.push_back( COO { res_col, o.o} );

                                    for ( auto &oi : currentOwnsO ) { if ( oi.o != o.o ) {

                                    double dd{ ( o.pnt - oi.o->pos ) . length() - oi.o->rad };

                                    if ( dd < oi.o->rad  ) {

                                        dd /= oi.o->rad;

                                       ray_res[ray_res.size()-1].colour =
                                                ( ( ray_res[ray_res.size()-1].colour ) * ( dd ) )
                                                + ( ( oi.o->col + ray_res[ray_res.size()-1].colour ) ) * ( 1 - dd );

                                        ray_res[ray_res.size()-1].colour.x = (ray_res[ray_res.size()-1].colour.x < 0 ) ? 0 : ( ray_res[ray_res.size()-1].colour.x > 255 ) ? 255 : ray_res[ray_res.size()-1].colour.x;
                                        ray_res[ray_res.size()-1].colour.y = (ray_res[ray_res.size()-1].colour.y < 0 ) ? 0 : ( ray_res[ray_res.size()-1].colour.y > 255 ) ? 255 : ray_res[ray_res.size()-1].colour.y;
                                        ray_res[ray_res.size()-1].colour.z = (ray_res[ray_res.size()-1].colour.z < 0 ) ? 0 : ( ray_res[ray_res.size()-1].colour.z > 255 ) ? 255 : ray_res[ray_res.size()-1].colour.z;

                                   } } } } } break; }

                  case oType::plane : {

                               for ( auto *o : bSphere->owns ) { double d1,d2;


                                 if ( o->intersect( orig, dir, d1, d2, res_col ) )
                                    { is_pln0 = true;
                                      double opx;

                                      opx  = ( d1 <= 20 ) ? d1 / 20.0 : ( d1 >= 227 ) ? ( 247 - d1 ) / 20.0 : 1.0 ;
                                      opx *= ( d2 <= 15 ) ? d2 / 15.0 : ( d2 >= 154 ) ? ( 169 - d2 ) / 15.0 : 1.0 ;

                                      ray_res.push_back( COO { res_col, o, opx } );

                                    } }

                           break; }

                  case oType::line : {

                                 for ( auto *o : bSphere->owns ) { double d1,d2;

                                   if ( o->intersect( orig, dir, d1, d2, res_col ) )

                                     ray_res.push_back( COO { res_col, o, d2 } );

                                  } for ( auto b : bSphere->bSpheres ) {

                                     it_Objects it; it( ray_res, b, dir, is1opq );
                                 }

                          break; }

                  case oType::point : { double d1, d2;

                        if  ( is1opq < 1.9 ) {
                            if ( !bSphere->owns.empty() ) {

                                 for ( auto *o : bSphere->owns ) {

                                   if ( o->intersect( orig, dir, d1, d2, res_col ) )

                                       pre_cdo.push_back( CDO { res_col, d1, d2 } );  }

                                 struct { bool operator()( CDO a, CDO b ) const { return a.dist < b.dist; } } sortCDO;
                                 std::sort( pre_cdo.begin(), pre_cdo.end(), sortCDO );

                                 for ( auto p : pre_cdo ) { if ( p.dist != INFINITY && is1opq < 1.9 )

                                     { ray_res.push_back( COO { p.colour, transObj, p.opq } ); is1opq += p.opq; }

                                 }  } if ( !bSphere->bSpheres.empty() && is1opq < 1.9 ) {

                                 struct bspc{ class bSphere* b0; double d0{INFINITY};
                                              bspc( class bSphere* bb0, double dd0 ) : b0{bb0}, d0{dd0} {} };

                                 std::vector<bspc> bvc;

                                 for ( auto b : bSphere->bSpheres )
                                 bvc.push_back( bspc { b, ( orig - b->pos ) . norm() } );

                                 struct { bool operator()( bspc a, bspc b ) const { return a.d0 < b.d0; } } sortbspc;
                                 std::sort( bvc.begin(), bvc.end(), sortbspc );

                                 for ( auto b : bvc ) { if ( b.b0->intersect( orig, dir, d1, d2, res_col) )

                                     {
                                         if ( b.b0->isvis ) ray_res.push_back( COO { Vec3d{0,33,99}, b.b0, 1 } ); //is1opq += b.b0->opq;
                                         it_Objects it; it( ray_res, b.b0, dir, is1opq );
                                         if ( b.b0->isvis ) ray_res.push_back( COO { Vec3d{0,33,99}, b.b0, 1 } ); //is1opq += b.b0->opq;

                                     }
                                 }

                            } }

                          break; }

                  case oType::order2 : {

                                for ( auto o : bSphere->owns2 ) { double d1{}, d2{}, o1{}, o2{}; Vec3d res_col2;

                                     if ( o->intersect( orig, dir, d1, o1, res_col, d2, o2, res_col2 ) )

                                     { pre_cdo.push_back( CDO { res_col, d1, o1 } );
                                       pre_cdo.push_back( CDO { res_col2, d2, o2 } ); }

                                 }

                                 struct { bool operator()( CDO a, CDO b ) const { return a.dist < b.dist; } } sortCDO;
                                 std::sort( pre_cdo.begin(), pre_cdo.end(), sortCDO );

                                 for ( auto p : pre_cdo ) { if ( p.dist != INFINITY )

                                     ray_res.push_back( COO { p.colour, transObj, p.opq } ); //std::clamp()

                                 }

                                 for ( auto b : bSphere->bSpheres ) {

                                     it_Objects it; it( ray_res, b, dir, is1opq );
                                 }


                          break; }

                        } // switch

                        } } itObjects; itObjects.operator()( ray_res, b.s, dir, is1opq );

                        if ( b.s->isvis ) ray_res.push_back( COO { b.s->col, b.s } );

                        newLen = INFINITY; }

                   } while ( newLen < prevLen ); if ( is1opq > 0.96 ) break;

                }

            }

        } itBVols;

       itBVols( ray_res, &currentBVols, dir, is1opq );

       Vec3d res_col;
       worldSphere.intersect( from, dir, pX, pY, res_col );
       ray_res.push_back( COO { res_col, &worldSphere } );

       if ( ray_res.size() == 0 ) { std::cout << "[ unexpected ] - /n"; dtX[p+0] = 0; dtX[p+1] = 0; dtX[p+2] = 0;  continue; }

       if ( ray_res.size() == 1 )
       { dtX[p+0] = ray_res[0].colour.x; dtX[p+1] = ray_res[0].colour.y; dtX[p+2] = ray_res[0].colour.z; continue; }

       if ( ray_res.size() >= 2 ) {

         if ( cursor ) pntd = ray_res[1].opnt;

         res_col = {0, 0, 0}; Vec3d t_col {0, 0, 0};

         double prevOpq{1.0}, tempOpq{0};

         for ( auto r : ray_res ) { prevOpq = 1 - tempOpq;

             if ( r.opq == 1.0 ) { tempOpq += r.opnt->opq;

             t_col = r.colour * r.opnt->opq * prevOpq;

             } else { t_col = r.colour * prevOpq * r.opq * r.opnt->opq;

                if ( r.opnt != &ss0 ) tempOpq += prevOpq * r.opq * r.opnt->opq; }

             res_col = res_col + t_col;

             if ( tempOpq > 0.99 ) break; }

         res_col.x = ( res_col.x > 255 ) ? 255 : res_col.x;
         res_col.y = ( res_col.y > 255 ) ? 255 : res_col.y;
         res_col.z = ( res_col.z > 255 ) ? 255 : res_col.z;

         dtX[p+0] = res_col.x; dtX[p+1] = res_col.y; dtX[p+2] = res_col.z; }

      }

   };

    int yL0{0};

    while ( yL0 < 780 ) { std::thread *th; int tcnt = ( 780 - yL0 > 12 ) ? num_threads : 780 - yL0;

        for ( int i = 0; i < tcnt; i++ )
        { th = new std::thread { processLine, yL0+i }; tcoll.push_back(th); }

        for ( int i = 0; i < tcnt; i++ )
        { tcoll[i]->join(); delete tcoll[i]; } tcoll.clear();

     yL0 += num_threads;

    }

    degtorad = qDegreesToRadians(phi);

    if ( bsp0.vis ) for ( auto t5 : tetra5 ) { t5->rotate( bsp00.pos, cos( degtorad ), sin( degtorad ), 0);
                                               t5->rotate( bsp00.pos, cos( degtorad ), sin( degtorad ), 1);
                                               t5->rotate( bsp00.pos, cos( degtorad ), sin( degtorad ), 2); }

    if ( bsp1.vis ) rotateSph();

    if ( is_pln0 )
    { planeTxt( false ); blur();

        rcnt++; if ( rcnt == 128 ) { rcnt = 0; ris = !ris; }
        gcnt++; if ( gcnt == 128 ) { gcnt = 0; gis = !gis; }
        bcnt++; if ( bcnt == 128 ) { bcnt = 0; bis = !bis; }
    }

    qDegRadX = qDegreesToRadians(h1x.angx); qDegRadY = qDegreesToRadians(h1x.angx); qDegRadZ = qDegreesToRadians(h1x.angx);

    h1x.rmx.m[0][0] = qCos(qDegRadX) * qCos(qDegRadY); h1x.rmx.m[0][1] = ( qCos(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) - ( qSin(qDegRadX) * qCos(qDegRadZ) ); h1x.rmx.m[0][2] = ( qCos(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) + ( qSin(qDegRadX) * qSin(qDegRadZ) );
    h1x.rmx.m[1][0] = qSin(qDegRadX) * qCos(qDegRadY); h1x.rmx.m[1][1] = ( qSin(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) + ( qCos(qDegRadX) * qCos(qDegRadZ) ); h1x.rmx.m[1][2] = ( qSin(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) - ( qCos(qDegRadX) * qSin(qDegRadZ) );
    h1x.rmx.m[2][0] = -1.0 * qSin(qDegRadY); h1x.rmx.m[2][1] = qCos(qDegRadY) * qSin(qDegRadZ); h1x.rmx.m[2][2] = qCos(qDegRadY) * qCos(qDegRadZ);

    h1x.angx += 0.078;

    qDegRadX = qDegreesToRadians(h2x.angx); qDegRadY = qDegreesToRadians(h2x.angx); qDegRadZ = qDegreesToRadians(h2x.angx);

    h2x.rmx.m[0][0] = qCos(qDegRadX) * qCos(qDegRadY); h2x.rmx.m[0][1] = ( qCos(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) - ( qSin(qDegRadX) * qCos(qDegRadZ) ); h2x.rmx.m[0][2] = ( qCos(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) + ( qSin(qDegRadX) * qSin(qDegRadZ) );
    h2x.rmx.m[1][0] = qSin(qDegRadX) * qCos(qDegRadY); h2x.rmx.m[1][1] = ( qSin(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) + ( qCos(qDegRadX) * qCos(qDegRadZ) ); h2x.rmx.m[1][2] = ( qSin(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) - ( qCos(qDegRadX) * qSin(qDegRadZ) );
    h2x.rmx.m[2][0] = -1.0 * qSin(qDegRadY); h2x.rmx.m[2][1] = qCos(qDegRadY) * qSin(qDegRadZ); h2x.rmx.m[2][2] = qCos(qDegRadY) * qCos(qDegRadZ);

    h2x.angx += 0.078;

    qDegRadX = qDegreesToRadians(h2y.angx); qDegRadY = qDegreesToRadians(h2y.angx); qDegRadZ = qDegreesToRadians(h2y.angx);

    h2y.rmx.m[0][0] = qCos(qDegRadX) * qCos(qDegRadY); h2y.rmx.m[0][1] = ( qCos(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) - ( qSin(qDegRadX) * qCos(qDegRadZ) ); h2y.rmx.m[0][2] = ( qCos(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) + ( qSin(qDegRadX) * qSin(qDegRadZ) );
    h2y.rmx.m[1][0] = qSin(qDegRadX) * qCos(qDegRadY); h2y.rmx.m[1][1] = ( qSin(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) + ( qCos(qDegRadX) * qCos(qDegRadZ) ); h2y.rmx.m[1][2] = ( qSin(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) - ( qCos(qDegRadX) * qSin(qDegRadZ) );
    h2y.rmx.m[2][0] = -1.0 * qSin(qDegRadY); h2y.rmx.m[2][1] = qCos(qDegRadY) * qSin(qDegRadZ); h2y.rmx.m[2][2] = qCos(qDegRadY) * qCos(qDegRadZ);

    h2y.angx++;


    l20.scw = 8 + double( control->ui->sld_sc->value() - 60 ) / 10.0;
    l21.scw = 8 + double( control->ui->sld_sc->value() - 60 ) / 10.0;
    l22.scw = 8 + double( control->ui->sld_sc->value() - 60 ) / 10.0;
    l23.scw = 8 + double( control->ui->sld_sc->value() - 60 ) / 10.0;

    h1x.scw = 1 - double( control->ui->sld_sc->value() - 79 ) / 150.0;
    h2x.scw = 1 - double( control->ui->sld_sc->value() - 79 ) / 150.0;

    // DOTS

    if ( bspT.vis ) {
    Vec3d pos{}, ppp{};
    Matrix44d a; uint cc{}, nn{};

    for ( auto p0 : pntarr ) { if ( cc == 39 ) { nn++; pos = {1,1,1}; ppp = {}; cc = 0; } cc++;

        double qDegRadX = qDegreesToRadians(p0->deg);
        double qDegRadY = qDegreesToRadians(p0->deg);
        double qDegRadZ = qDegreesToRadians(p0->deg);

        double sc = 1.15 - std::abs( 15.0 - p0->deg );
        if ( !p0->yn ) { if ( p0->deg > 14.0 ) p0->deg -= 0.25 * sc; else { p0->yn = true; p0->deg = 14.0; } }
        else { if ( p0->deg < 16.0 ) p0->deg += 0.25 * sc; else { p0->yn = false; p0->deg = 16.0; } }

        a.m[0][0] = qCos(qDegRadX) * qCos(qDegRadY); a.m[0][1] = ( qCos(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) - ( qSin(qDegRadX) * qCos(qDegRadZ) ); a.m[0][2] = ( qCos(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) + ( qSin(qDegRadX) * qSin(qDegRadZ) );
        a.m[1][0] = qSin(qDegRadX) * qCos(qDegRadY);  a.m[1][1] = ( qSin(qDegRadX) * qSin(qDegRadY) * qSin(qDegRadZ) ) + ( qCos(-qDegRadX) * qCos(qDegRadZ) ); a.m[1][2] = ( qSin(qDegRadX) * qSin(qDegRadY) * qCos(qDegRadZ) ) - ( qCos(qDegRadX) * qSin(qDegRadZ) );
        a.m[2][0] = -1.0 * qSin(qDegRadY);            a.m[2][1] = qCos(qDegRadY) * qSin(qDegRadZ);                                                            a.m[2][2] = qCos(qDegRadY) * qCos(qDegRadZ);

        p0->pos = ppp + bspT.pos;

        pos = pos * 1.087;
        a.multDirMatrix( pos, ppp );
        pos = ppp;

        if ( p0->owner && ( p0->pos - p0->owner->pos ).length() > p0->owner->rad ) { //p0->owner->isvis = 1;

        p0->owner->rad = ( p0->pos - p0->owner->pos ).length() + 15;
        p0->owner->rad2 = p0->owner->rad * p0->owner->rad;

        if ( p0->owner->owner->rad < p0->owner->rad ) { p0->owner->owner->rad2 = p0->owner->rad2; }

        }

        switch ( nn ) {

        case 0 : break;
        case 1 : ppp.x *= -1; break;
        case 2 : ppp.y *= -1; break;
        case 3 : ppp.z *= -1; break;
        case 4 : ppp.x *= -1; ppp.z *= -1; break;
        case 5 : ppp.y *= -1; ppp.z *= -1; break;
        case 6 : ppp.x *= -1; ppp.y *= -1; break;
        case 7 : ppp.x *= -1; ppp.y *= -1; ppp.z *= -1; break;

       }
    }



    Vec3d bpos { bspT.pos }; uint icnt{16};

    while ( icnt != 2 ) {

        createInbSpheres( bpos );

    switch ( icnt ) {

    case 16 : bpos = bspT.pos + 13; break;
    case 15 : bpos = bspT.pos - 13; break;
    case 14 : bpos = bspT.pos + 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; break;
    case 13 : bpos = bspT.pos - 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; break;

    case 12 : bpos = bspT.pos - 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 11 : bpos = bspT.pos + 13; bpos.z = ( ( bpos.z - bspT.pos.z ) * -1 ) + bspT.pos.z; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 10 : bpos = bspT.pos - 13; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;
    case 9  : bpos = bspT.pos + 13; bpos.x = ( ( bpos.x - bspT.pos.x ) * -1 ) + bspT.pos.x; break;

    case 8 : bpos = bspT.pos; bpos.z += 23; break;
    case 7 : bpos = bspT.pos; bpos.z -= 23; break;
    case 6 : bpos = bspT.pos; bpos.y += 23; break;
    case 5 : bpos = bspT.pos; bpos.y -= 23; break;
    case 4 : bpos = bspT.pos; bpos.x += 23; break;
    case 3 : bpos = bspT.pos; bpos.x -= 23; break;
    }

    icnt--; } createInbSpheres( bpos );

    for ( auto p0 : pntarr ) { bool isin { false };

        for ( auto b0 : bspT.bSpheres ) {

            if ( ( p0->pos - b0->pos ) . length() < b0->rad - 1 ) {

                for ( auto b1 : b0->bSpheres ) {

                    if ( ( p0->pos - b1->pos ) . length() < b1->rad - 3.1 ) {

                        if ( !p0->owner ) { b1->owns.push_back(p0); isin = true; p0->owner = b1; }

                    } if ( isin ) break;

                } if ( !isin && !p0->owner ) { std::cout << "missed <- \n"; b0->owns.push_back(p0);
                       p0->owner = b0; p0->col = { 255, 255, 255 }; isin = true; }

            }
        }
    }

    for ( auto b0 : bspT.bSpheres ) {

        for ( std::vector<bSphere*>::iterator it = b0->bSpheres.begin(); it != b0->bSpheres.end();) {

            if ( it.operator*()->owns.empty() ) { b0->bSpheres.erase( it ); }
            else {  Vec3d newpos{}; double ocnt{};
                    for ( auto o : it.operator*()->owns ) {
                    newpos = newpos + ( o->pos - it.operator*()->pos );
                    ocnt++;
                    } newpos = { newpos.x / ocnt, newpos.y / ocnt, newpos.z / ocnt };
                    it.operator*()->pos = newpos + it.operator*()->pos;
                    double res{};
                    for ( auto o : it.operator*()->owns ) {
                    if ( res < ( o->pos - it.operator*()->pos ) . length() + 1  )
                       { res = ( o->pos - it.operator*()->pos ) . length() + 1; }
                    } it.operator*()->rad = res; it.operator*()->rad2 = res * res;
                    if ( res > 10.5 ) std::cout << res << "\n";
                    it++; }

        } if ( !b0->owns.empty() && b0->owns.size() > 1 ) { std::cout << "not empty - \n"; }

    }  for ( std::vector<bSphere*>::iterator it = bspT.bSpheres.begin(); it != bspT.bSpheres.end();)
           { if ( it.operator*()->bSpheres.empty() ) bspT.bSpheres.erase(it);
             else { Vec3d newpos{  }; double ocnt{};
                    for ( auto b : it.operator*()->bSpheres ) {
                    newpos = newpos + ( b->pos - it.operator*()->pos ); ocnt++;
                    } newpos = { newpos.x / ocnt, newpos.y / ocnt, newpos.z / ocnt };
                    it.operator*()->pos = newpos + it.operator*()->pos;
                    double res{};
                    for ( auto b : it.operator*()->bSpheres ) {
                        if ( res < ( b->pos - it.operator*()->pos ) . length() + b->rad + 1 )
                           { res = ( b->pos - it.operator*()->pos ) . length() + b->rad + 1; }
                    } it.operator*()->rad = res; it.operator*()->rad2 = res * res;
                    it++; } }

    }

    looptime = timer.elapsed();

}

void MainWindow::endF()
{

    QString s, lx, ly, lz;

    if ( control->changed ) { cursis = control->cursor;

     double rx, ry, rz;



     rx = lx.toDouble();
     ry = ly.toDouble();
     rz = lz.toDouble();

          qDegRadX = qDegreesToRadians(rx); qDegRadY = qDegreesToRadians(ry); qDegRadZ = qDegreesToRadians(rz);

          std::cout << rx << " " << ry << " " << rz << '\n';

           Vec3d aa, bb, cc, a, b, c;
           for ( int i = 16; i <=19; i++ ) {

           }


        control->changed = false;
    }

    info->ui->txtMsgX->setText( "Tx: " + s.setNum( to.x ) );
    info->ui->txtMsgY->setText( "Ty: " + s.setNum( to.y ) );
    info->ui->txtMsgZ->setText( "Tz: " + s.setNum( to.z ) );
    info->ui->txtMsgT->setText( "T: " + s.setNum( looptime, 'f', 6 ) );

    worldSphere.col.x = 88 * control->ui->sldb->value() / 100;
    worldSphere.col.y = 66 * control->ui->sldg->value() / 100;
    worldSphere.col.z = 77 * control->ui->sldr->value() / 100;

   pix0->convertFromImage( *img0, 0 );
   qGpixItem0->setPixmap( *pix0 );

    resc = control->sld00_v;
    info->ui->txtMsgA->setText( "A: " + s.setNum( resc, 'f', 6 ) );

    info->ui->txtMsgOx->setText( "Ox: " + s.setNum( from.x, 'f', 3 ) );
    info->ui->txtMsgOy->setText( "Oy: " + s.setNum( from.y, 'f', 3 ) );
    info->ui->txtMsgOz->setText( "Oz: " + s.setNum( from.z, 'f', 3 ) );
    scwvec = control->ui->sld00->value() * resc;
    info->ui->txtMsgS->setText( "S: " + s.setNum( timer.elapsed(), 'f', 6 ) );

}

void MainWindow::rotateSph()
{
    ( sb0 ) ? sc0 -= 0.001 : sc0 += 0.001; if ( sc0 < 0.71 || sc0 > 1.01 ) sb0 = !sb0;
    ( sb1 ) ? sc1 -= 0.002 : sc1 += 0.002; if ( sc1 < 0.63 || sc1 > 1.09 ) sb1 = !sb1;
    ( sb2 ) ? sc2 -= 0.003 : sc2 += 0.003; if ( sc2 < 0.56 || sc2 > 1.03 ) sb2 = !sb2;
    ( sb3 ) ? sc3 -= 0.004 : sc3 += 0.004; if ( sc3 < 0.64 || sc3 > 1.06 ) sb3 = !sb3;

vsx.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc0)), qSin(qDegreesToRadians(3.99*sc0)), 2);
vsc.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99)), qSin(qDegreesToRadians(3.99)), 0); vsc.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc0)), qSin(qDegreesToRadians(3.99*sc0)), 2);
vsy.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc1)), qSin(qDegreesToRadians(3.99*sc1)), 1);
vsb.rotate(bsp1.pos, qCos(qDegreesToRadians(1.99)), qSin(qDegreesToRadians(1.99)), 0); vsb.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc1)), qSin(qDegreesToRadians(3.99*sc1)), 1);
vsw.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc3)), qSin(qDegreesToRadians(3.99*sc2)), 2);
vsn.rotate(bsp1.pos, qCos(qDegreesToRadians(1.99)), qSin(qDegreesToRadians(1.99)), 1); vsn.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc3)), qSin(qDegreesToRadians(3.99*sc2)), 2);

vsq.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc2)), qSin(qDegreesToRadians(3.99*sc2)), 0);
vsi.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99)), qSin(qDegreesToRadians(3.99)), 1); vsi.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc2)), qSin(qDegreesToRadians(3.99*sc2)), 2);
vsr.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc1)), qSin(qDegreesToRadians(3.99*sc1)), 2);
vsg.rotate(bsp1.pos, qCos(qDegreesToRadians(1.69)), qSin(qDegreesToRadians(1.69)), 2); vsg.rotate(bsp1.pos, qCos(qDegreesToRadians(2.39*sc1)), qSin(qDegreesToRadians(2.39*sc1)), 1);
vsf.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc2)), qSin(qDegreesToRadians(3.99*sc2)), 1);
vss.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99)), qSin(qDegreesToRadians(3.99)), 2); vss.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc2)), qSin(qDegreesToRadians(3.99*sc2)), 1);

vsd.rotate(bsp1.pos, qCos(qDegreesToRadians(6.69*sc3)), qSin(qDegreesToRadians(6.69*sc3)), 1);
vsa.rotate(bsp1.pos, qCos(qDegreesToRadians(1.99)), qSin(qDegreesToRadians(1.99)), 2); vsa.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc3)), qSin(qDegreesToRadians(3.99*sc3)), 0);
vsz.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc0)), qSin(qDegreesToRadians(3.99*sc0)), 0);
vse.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99)), qSin(qDegreesToRadians(3.99)), 1); vse.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc0)), qSin(qDegreesToRadians(3.99*sc0)), 2);
vsm.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc1)), qSin(qDegreesToRadians(3.99*sc1)), 2);
vsv.rotate(bsp1.pos, qCos(qDegreesToRadians(1.99)), qSin(qDegreesToRadians(1.99)), 0); vsv.rotate(bsp1.pos, qCos(qDegreesToRadians(3.99*sc1)), qSin(qDegreesToRadians(3.99*sc1)), 2);

ss1b.pos = ( ( vsb - bsp1.pos ) * sc3 ) + bsp1.pos; ss1c.pos = ( ( vsc - bsp1.pos ) * sc1 ) + bsp1.pos;
ss1g.pos = ( ( vsg - bsp1.pos ) * sc2 ) + bsp1.pos; ss1e.pos = ( ( vse - bsp1.pos ) * sc0 ) + bsp1.pos;
ss1i.pos = ( ( vsi - bsp1.pos ) * sc1 ) + bsp1.pos; ss1n.pos = ( ( vsn - bsp1.pos ) * sc0 ) + bsp1.pos;

ss1q.pos = ( ( vsq - bsp1.pos ) * sc3 ) + bsp1.pos; ss1w.pos = ( ( vsw - bsp1.pos ) * sc1 ) + bsp1.pos;
ss1r.pos = ( ( vsr - bsp1.pos ) * sc2 ) + bsp1.pos; ss1y.pos = ( ( vsy - bsp1.pos ) * sc0 ) + bsp1.pos;
ss1f.pos = ( ( vsf - bsp1.pos ) * sc2 ) + bsp1.pos; ss1s.pos = ( ( vss - bsp1.pos ) * sc1 ) + bsp1.pos;

ss1d.pos = ( ( vsd - bsp1.pos ) * sc3 ) + bsp1.pos; ss1a.pos = ( ( vsa - bsp1.pos ) * sc1 ) + bsp1.pos;
ss1z.pos = ( ( vsz - bsp1.pos ) * sc2 ) + bsp1.pos; ss1x.pos = ( ( vsx - bsp1.pos ) * sc0 ) + bsp1.pos;
ss1m.pos = ( ( vsm - bsp1.pos ) * sc2 ) + bsp1.pos; ss1v.pos = ( ( vsv - bsp1.pos ) * sc3 ) + bsp1.pos;
}

void MainWindow::planeTxt(bool is)
{

for ( int x = 1; x <= 494; x++ ) { for ( int y = 1; y <= 338; y++ ) {

    int p = ( ( y * 494 ) + x ) * 3, pt;

    if ( !bis ) dtO[p+0] -= 1; else dtO[p+0] += 1;
    if ( !gis ) dtO[p+1] += 1; else dtO[p+1] -= 1;
    if ( !ris ) dtO[p+2] += 1; else dtO[p+2] -= 1;

    Vec3d res = { (double)dtO[p+0], (double)dtO[p+1], (double)dtO[p+2] };

    pt = ( ( ( y - 0 ) * 494 ) + ( x - 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y - 1 ) * 494 ) + ( x - 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y + 1 ) * 494 ) + ( x - 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y - 1 ) * 494 ) + ( x - 0 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y + 1 ) * 494 ) + ( x - 0 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y + 0 ) * 494 ) + ( x + 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y - 1 ) * 494 ) + ( x + 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };
    pt = ( ( ( y + 1 ) * 494 ) + ( x + 1 ) ) * 3;
    res = res + Vec3d{ (double)dtO[pt+0], (double)dtO[pt+1], (double)dtO[pt+2] };

    res = res * 0.111111111111;

    dtZ[p+0] = res.x; dtZ[p+1] = res.y; dtZ[p+2] = res.z;

    } }

    if ( is ) { for ( int x = 1; x <= 494; x++ ) { for ( int y = 1; y <= 338; y++ ) {

        int p = ( ( y * 494 ) + x ) * 3;

        dtO[p+0] = dtZ[p+0]; dtO[p+1] = dtZ[p+1]; dtO[p+2] = dtZ[p+2];

  } } } else { for ( int x = 1; x <= 494; x++ ) { for ( int y = 1; y <= 338; y++ ) {

            int p = ( ( y * 494 ) + x ) * 3, pt; double p_x, p_y, res;

            Vec3d a00,b00,c00,a01,b01,c01,a02,b02,c02;

            pt = ( ( ( y - 1 ) * 494 ) + ( x - 1 ) ) * 3;
            a00 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y - 0 ) * 494 ) + ( x - 1 ) ) * 3;
            b00 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y + 1 ) * 494 ) + ( x - 1 ) ) * 3;
            c00 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };

            pt = ( ( ( y - 1 ) * 494 ) + ( x - 0 ) ) * 3;
            a01 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y - 0 ) * 494 ) + ( x - 0 ) ) * 3;
            b01 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y + 1 ) * 494 ) + ( x - 0 ) ) * 3;
            c01 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };

            pt = ( ( ( y - 1 ) * 494 ) + ( x + 1 ) ) * 3;
            b02 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y + 0 ) * 494 ) + ( x + 1 ) ) * 3;
            a02 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };
            pt = ( ( ( y + 1 ) * 494 ) + ( x + 1 ) ) * 3;
            c02 = Vec3d{ (double)dtZ[pt+0], (double)dtZ[pt+1], (double)dtZ[pt+2] };


            p_x = (ksx[0][0] * a00.x) + (ksx[0][1] * a01.x) + (ksx[0][2] * a02.x) +
                  (ksx[1][0] * b00.x) + (ksx[1][1] * b01.x) + (ksx[1][2] * b02.x) +
                  (ksx[2][0] * c00.x) + (ksx[2][1] * c01.x) + (ksx[2][2] * c02.x) ;

            p_y = (ksy[0][0] * a00.x) + (ksy[0][1] * a01.x) + (ksy[0][2] * a02.x) +
                  (ksy[1][0] * b00.x) + (ksy[1][1] * b01.x) + (ksy[1][2] * b02.x) +
                  (ksy[2][0] * c00.x) + (ksy[2][1] * c01.x) + (ksy[2][2] * c02.x) ;

            res = sqrt ( ( p_x * p_x ) + ( p_y * p_y ) );

            res += ( dtZ[p+0] & 0b00111111 ) + ( dtZ[p+1] & 0b00000111 ) + ( dtZ[p+2] & 0b10011001 );

            res = ( res < 0 ) ? 0 : ( res > 255 ) ? 255 : res;

            dtJ[p+0] =  res;

            p_x = (ksx[0][0] * a00.y) + (ksx[0][1] * a01.y) + (ksx[0][2] * a02.y) +
                  (ksx[1][0] * b00.y) + (ksx[1][1] * b01.y) + (ksx[1][2] * b02.y) +
                  (ksx[2][0] * c00.y) + (ksx[2][1] * c01.y) + (ksx[2][2] * c02.y) ;

            p_y = (ksy[0][0] * a00.y) + (ksy[0][1] * a01.y) + (ksy[0][2] * a02.y) +
                  (ksy[1][0] * b00.y) + (ksy[1][1] * b01.y) + (ksy[1][2] * b02.y) +
                  (ksy[2][0] * c00.y) + (ksy[2][1] * c01.y) + (ksy[2][2] * c02.y) ;

            res = sqrt ( ( p_x * p_x ) + ( p_y * p_y ) );

            res += + ( dtZ[p+1] & 0b00111111 ) + ( dtZ[p+2] & 0b10000001 ) + ( dtZ[p+0] & 0b00111000 );

            res = ( res < 0 ) ? 0 : ( res > 255 ) ? 255 : res;

            dtJ[p+1] = res ;

            p_x = (ksx[0][0] * a00.z) + (ksx[0][1] * a01.z) + (ksx[0][2] * a02.z) +
                  (ksx[1][0] * b00.z) + (ksx[1][1] * b01.z) + (ksx[1][2] * b02.z) +
                  (ksx[2][0] * c00.z) + (ksx[2][1] * c01.z) + (ksx[2][2] * c02.z) ;

            p_y = (ksy[0][0] * a00.z) + (ksy[0][1] * a01.z) + (ksy[0][2] * a02.z) +
                  (ksy[1][0] * b00.z) + (ksy[1][1] * b01.z) + (ksy[1][2] * b02.z) +
                  (ksy[2][0] * c00.z) + (ksy[2][1] * c01.z) + (ksy[2][2] * c02.z) ;

            res = sqrt ( ( p_x * p_x ) + ( p_y * p_y ) );

            res += ( dtZ[p+2] & 0b00111111 ) + ( dtZ[p+1] & 0b10000001 ) + ( dtZ[p+0] & 0b00101010 );

            res = ( res < 0 ) ? 0 : ( res > 255 ) ? 255 : res;

            dtJ[p+2] = res ;
    } } }

}

void MainWindow::blur()
{

  for ( int x = 1; x <= 494; x++ ) { for ( int y = 1; y <= 338; y++ ) {

                int p = ( ( y * 494 ) + x ) * 3, pt;

                Vec3d res = { (double)dtJ[p+0], (double)dtJ[p+1], (double)dtJ[p+2] };

                pt = ( ( ( y - 0 ) * 494 ) + ( x - 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y - 1 ) * 494 ) + ( x - 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y + 1 ) * 494 ) + ( x - 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y - 1 ) * 494 ) + ( x - 0 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y + 1 ) * 494 ) + ( x - 0 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y + 0 ) * 494 ) + ( x + 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y - 1 ) * 494 ) + ( x + 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };
                pt = ( ( ( y + 1 ) * 494 ) + ( x + 1 ) ) * 3;
                res = res + Vec3d{ (double)dtJ[pt+0], (double)dtJ[pt+1], (double)dtJ[pt+2] };

                res = res * 0.111111111;

                dtJ[p+0] = res.x; dtJ[p+1] = res.y; dtJ[p+2] = res.z;

                } }
}

void MainWindow::createInbSpheres( Vec3d bpos )
{
    bSphere* bspT9 = new bSphere( bpos, ( ( bpos + 5 ) - bpos ) . length() + 7, 0, &bspT, oType::point ); //rbSph.push_back( bspT9 );

    bSphere* bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );

    bspT0 = new bSphere( bpos-5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0 = new bSphere( bpos+5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0 = new bSphere( bpos-5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.z = ( ( bspT0->pos.z - bpos.z ) * -1.0 ) + bpos.z;
    bspT0 = new bSphere( bpos+5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.z = ( ( bspT0->pos.z - bpos.z ) * -1.0 ) + bpos.z;

    bspT0 = new bSphere( bpos-5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.z = ( ( bspT0->pos.z - bpos.z ) * -1.0 ) + bpos.z;
    bspT0->pos.x = ( ( bspT0->pos.x - bpos.x ) * -1.0 ) + bpos.x;
    bspT0 = new bSphere( bpos+5, 12.3, 0, bspT9, oType::point );// rbSph.push_back( bspT0 );
    bspT0->pos.z = ( ( bspT0->pos.z - bpos.z ) * -1.0 ) + bpos.z;
    bspT0->pos.x = ( ( bspT0->pos.x - bpos.x ) * -1.0 ) + bpos.x;
    bspT0 = new bSphere( bpos-5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.x = ( ( bspT0->pos.x - bpos.x ) * -1.0 ) + bpos.x;
    bspT0 = new bSphere( bpos+5, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.x = ( ( bspT0->pos.x - bpos.x ) * -1.0 ) + bpos.x;

    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.z += 8.5;
    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.z -= 8.5;
    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point );// rbSph.push_back( bspT0 );
    bspT0->pos.y += 8.5;
    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point );// rbSph.push_back( bspT0 );
    bspT0->pos.y -= 8.5;
    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point );// rbSph.push_back( bspT0 );
    bspT0->pos.x += 8.5;
    bspT0 = new bSphere( bpos, 12.3, 0, bspT9, oType::point ); //rbSph.push_back( bspT0 );
    bspT0->pos.x -= 8.5;

}
