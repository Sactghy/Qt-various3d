#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <iostream>
#include <iomanip>
#include <vector>

#include <QtMath>

static double scwvec{1};
constexpr double kEpsilon = 1e-7;
enum class oType { none, sphere, triangle, plane, line, order2, point };

template<typename T>
class Vec3
{ public: T x, y, z;
    Vec3() : x(T(1)), y(T(1)), z(T(1)) {}
    Vec3(const T &xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}

    inline T norm() const { return x * x + y * y + z * z; }

    inline T length() const { return sqrt( x * x + y * y + z * z ); }

    Vec3<T>& normalize()
    { T len = norm();
     // if (len > 0) {
         T invLen = 1 / std::sqrt(len);
         x *= invLen, y *= invLen, z *= invLen;
      //}
        return *this; }

    inline T dot(const Vec3<T> &v) const
    { return x * v.x + y * v.y + z * v.z; }

    Vec3<T> cross(const Vec3<T> &v) const
    { return Vec3<T>( y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }

    Vec3<T> operator + (const Vec3<T> &v) const
        { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T> operator - (const Vec3<T> &v) const
        { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator * (const T &r) const
        { return Vec3<T>(x * r, y * r, z * r); }
    Vec3<T> operator * (const Vec3<T> &r) const
        { return Vec3<T>(x * r.x, y * r.y, z * r.z); }

    const T& operator [] (int i) const { return (&x)[i]; }
    T& operator [] (int i) { return (&x)[i]; }

    friend std::ostream& operator << (std::ostream &s, const Vec3<T> &v)
    { return s << "v " << v.x << ' ' << v.y << ' ' << v.z << '\n'; }

    void rotate(Vec3<T> ps, double cx, double sy, int tp)

         {   double tx, tz;

             if ( tp == 1 ) { x = x - ps.x; y = y - ps.y; z = z - ps.z; tx=(x*cx)-(z*sy); tz=(z*cx)+(x*sy);
                              z = tz; x = tx; x = x + ps.x; y = y + ps.y; z = z + ps.z; }

             if ( tp == 0 ) { x = x - ps.x; y = y - ps.y; z = z - ps.z; tx=(x*cx)-(y*sy); tz=(y*cx)+(x*sy);
                              y = tz; x = tx; x = x + ps.x; y = y + ps.y; z = z + ps.z; }

             if ( tp == 2 ) { x = x - ps.x; y = y - ps.y; z = z - ps.z; tx=(y*cx)-(z*sy); tz=(z*cx)+(y*sy);
                              z = tz; y = tx; x = x + ps.x; y = y + ps.y; z = z + ps.z; }
         }


    friend bool operator< (const Vec3& a1, const Vec3& a2)
    { return ( ( a1.x < a2.x ? true : false ) & ( a1.y < a2.y ? true : false ) & ( a1.z < a2.z ? true : false ) ); }
    friend bool operator> (const Vec3& a1, const Vec3& a2)
    { return ( ( a1.x > a2.x ? true : false ) & ( a1.y > a2.y ? true : false ) & ( a1.z > a2.z ? true : false ) ); }

    friend bool operator<= (const Vec3& a1, const Vec3& a2)
    { return ( ( a1.x <= a2.x ? true : false ) || ( a1.y <= a2.y ? true : false ) || ( a1.z <= a2.z ? true : false ) ); }
    friend bool operator>= (const Vec3& a1, const Vec3& a2)
    { return ( ( a1.x >= a2.x ? true : false ) || ( a1.y >= a2.y ? true : false ) || ( a1.z >= a2.z ? true : false ) ); }

    friend bool operator== (const Vec3& a1, const Vec3& a2)
    { return ( ( a1.x == a2.x ? true : false ) & ( a1.y == a2.y ? true : false ) & ( a1.z == a2.z ? true : false ) ); }

};

template<typename T>
class Vec2
{ public: T x, y;
    Vec2() : x( T(0) ), y( T(0) ) {}
    Vec2( const T &xx ) : x(xx), y(xx) {}
    Vec2( T xx, T yy ) : x(xx), y(yy) {}

    T norm() const { return x * x + y * y; }

    T length() const { return sqrt( x * x + y * y ); }

    Vec2<T>& normalize()
    { T len = norm();
      if (len > 0) { T invLen = 1 / sqrt(len);
         x *= invLen, y *= invLen; }
        return *this; }

    T dot(const Vec2<T> &v) const
    { return x * v.x + y * v.y; }

    Vec2<T> operator + (const Vec2<T> &v) const
        { return Vec2<T>( x + v.x, y + v.y ); }
    Vec2<T> operator - (const Vec2<T> &v) const
        { return Vec2<T>( x - v.x, y - v.y ); }
    Vec2<T> operator * (const T &r) const
        { return Vec2<T>( x * r, y * r ); }
    Vec2<T> operator * (const Vec2<T> &r) const
        { return Vec2<T>( x * r.x, y * r.y ); }

    const T& operator [] (int i) const { return (&x)[i]; }
    T& operator [] (int i) { return (&x)[i]; }

    friend std::ostream& operator << (std::ostream &s, const Vec2<T> &v)
    { return s << "v " << v.x << ' ' << v.y << '\n'; }

    void rotate(Vec2<T> ps, double cx, double sy, int tp)

         {   double tx, ty;

             if ( tp == 0 ) { x = x - ps.x; y = y - ps.y; tx=(x*cx)-(y*sy); ty=(x*sy)+(y*cx);
                              x = tx; y = ty; x = x + ps.x; y = y + ps.y; }

         }


};

template<typename T>
class Matrix44
{
public:
    T m[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

    Matrix44(){}

    Matrix44 (T a, T b, T c, T d, T e, T f, T g, T h, T i, T j, T k, T l, T m, T n, T o, T p)
    { m[0][0]=a; m[0][1]=b; m[0][2]=c; m[0][3]=d; m[1][0]=e; m[1][1]=f; m[1][2]=g; m[1][3]=h; m[2][0]=i; m[2][1]=j; m[2][2]=k; m[2][3]=l; m[3][0]=m; m[3][1]=n; m[3][2]=o; m[3][3]=p; }

    const T* operator [] (int i) const { return m[i]; }

    T* operator [] (int i) { return m[i]; }

    friend std::ostream& operator << (std::ostream &s, const Matrix44 &m)
    {   std::ios_base::fmtflags oldFlags = s.flags();
        int width = 12; s.precision(5); s.setf (std::ios_base::fixed);
        s << "(" << std::setw (width) << m[0][0] << " " << std::setw (width) << m[0][1] << " " << std::setw (width) << m[0][2] << " " << std::setw (width) << m[0][3] << "\n" <<
             " " << std::setw (width) << m[1][0] << " " << std::setw (width) << m[1][1] << " " << std::setw (width) << m[1][2] << " " << std::setw (width) << m[1][3] << "\n" <<
             " " << std::setw (width) << m[2][0] << " " << std::setw (width) << m[2][1] << " " << std::setw (width) << m[2][2] << " " << std::setw (width) << m[2][3] << "\n" <<
             " " << std::setw (width) << m[3][0] << " " << std::setw (width) << m[3][1] << " " << std::setw (width) << m[3][2] << " " << std::setw (width) << m[3][3] << ")\n";
        s.flags (oldFlags); return s; }

    Matrix44 operator * (const Matrix44& rhs) const
    {   Matrix44 mult;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                mult[i][j]=m[i][0]*rhs[0][j]+m[i][1]*rhs[1][j]+m[i][2]*rhs[2][j]+m[i][3]*rhs[3][j];  } }
        return mult;  }

    Matrix44 transposed() const
        { return Matrix44 (m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1], m[0][2], m[1][2], m[2][2], m[3][2], m[0][3], m[1][3], m[2][3], m[3][3]); }
    Matrix44& transpose ()
        { Matrix44 tmp (m[0][0], m[1][0], m[2][0], m[3][0], m[0][1], m[1][1], m[2][1], m[3][1], m[0][2], m[1][2], m[2][2], m[3][2], m[0][3], m[1][3], m[2][3], m[3][3]);
            *this = tmp; return *this; }

    void multVecMatrix(const Vec3<T> &src, Vec3<T> &dst) const
    {   T w;
        dst.x = src.x * m[0][0] + src.y * m[1][0] + src.z * m[2][0] + m[3][0];
        dst.y = src.x * m[0][1] + src.y * m[1][1] + src.z * m[2][1] + m[3][1];
        dst.z = src.x * m[0][2] + src.y * m[1][2] + src.z * m[2][2] + m[3][2];
        w = src.x * m[0][3] + src.y * m[1][3] + src.z * m[2][3] + m[3][3];
        if (w != 1 || w != 0) { dst.x /= w; dst.y /= w; dst.z /= w;  }   }

    void multDirMatrix(const Vec3<T> &src, Vec3<T> &dst) const
    {   T a, b, c;
        a = src[0] * m[0][0] + src[1] * m[1][0] + src[2] * m[2][0];
        b = src[0] * m[0][1] + src[1] * m[1][1] + src[2] * m[2][1];
        c = src[0] * m[0][2] + src[1] * m[1][2] + src[2] * m[2][2];
        dst.x = a; dst.y = b; dst.z = c; }

    Matrix44 inverse()
        {   int i, j, k; Matrix44 s; Matrix44 t (*this);
            // Forward elimination
            for (i = 0; i < 3 ; i++) { int pivot = i; T pivotsize = t[i][i];
                 if (pivotsize < 0) pivotsize = -pivotsize;
                 for (j = i + 1; j < 4; j++) { T tmp = t[j][i];
                     if (tmp < 0) tmp = -tmp;
                     if (tmp > pivotsize) { pivot = j; pivotsize = tmp; } }
                 if (pivotsize == 0) { // Cannot invert singular matrix
                     return Matrix44(); }
                 if (pivot != i) { for (j = 0; j < 4; j++) {
                        T tmp;
                        tmp = t[i][j];
                        t[i][j] = t[pivot][j];
                        t[pivot][j] = tmp;
                        tmp = s[i][j];
                        s[i][j] = s[pivot][j];
                        s[pivot][j] = tmp;  } }
                 for (j = i + 1; j < 4; j++) { T f = t[j][i] / t[i][i];
                        for (k = 0; k < 4; k++) { t[j][k] -= f * t[i][k]; s[j][k] -= f * s[i][k]; } }
                                    }
          // Backward substitution
            for (i = 3; i >= 0; --i) { T f;
                 if ((f = t[i][i]) == 0) { // Cannot invert singular matrix
                    return Matrix44(); }
                for (j = 0; j < 4; j++) { t[i][j] /= f; s[i][j] /= f; }
                for (j = 0; j < i; j++) { f = t[j][i];
                    for (k = 0; k < 4; k++) { t[j][k] -= f * t[i][k]; s[j][k] -= f * s[i][k]; } }
                            } return s; }

    const Matrix44<T>& invert() { *this = inverse(); return *this; }

};

typedef Matrix44<double> Matrix44d;
typedef Vec3<double> Vec3d;
typedef Vec2<double> Vec2d;
typedef Vec3<unsigned char> Vec3u;

class bSphere;

class Object
{
    public: Vec3d pos, col;
            bool vis{true}, isvis{}, blur{false};
            double opq{}, rad{};
            std::vector<bSphere*> bSpheres;
            oType otype{oType::none};
            Object* pair;

   virtual ~Object() {};

   virtual void scale (double scaler) = 0;

   virtual void rotate(Vec3d ps, double cx, double sy, int tp) = 0;

   virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) = 0;

};

class Object2
{
    public: Vec3d pos, col;
            bool vis{true}, isvis{}, blur{false};
            double opq{}, rad{};
            std::vector<bSphere*> bSpheres;
            oType otype{oType::none};
            Object* pair;

   virtual ~Object2() {};

   virtual void scale (double scaler) = 0;

   virtual void rotate( Vec3d ps, double cx, double sy, int tp ) = 0;

   virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &op1, Vec3d &res_col1, double &dist2, double &op2, Vec3d &res_col2) = 0;

};


class wSphere : public Object
{
    public: double rad2;
            bool type {0};
            Vec3d o_pos{};

    wSphere( Vec3d p = Vec3d(0,0,0), double r = 69, Vec3d c = Vec3d(0,0,0), double o = 1.0 )
    {
        pos = p; rad = r; rad2 = r * r; col = c; opq = o; isvis = true;
    }

    virtual void rotate(Vec3d, double, double, int ) override
    { }

   virtual void scale (double scaler) override
   {
     rad2 *= scaler * scaler;
     if ( o_pos.x || o_pos.y || o_pos.z) { pos = ( ( pos - o_pos ) * scaler ) + o_pos; }
   }

   inline virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
   {
       Vec3d L = pos - orig;
       double tca = L.dot(dir),
              thc = sqrt( rad2 - ( L.dot(L) - tca * tca )  );
       dist1 = tca - thc; dist2 = tca + thc;

       dist1 = ( dist1 < dist2 ) ?  dist2  :  dist1;

       Vec3d Phit = ( ( orig + dir * dist1 ) - pos ) . normalize() * scwvec;

       res_col = { col.x + ( Phit.x / M_PI ) * 100 ,
                   col.y - ( Phit.y / M_PI ) * 100,
                   col.z + ( Phit.z / M_PI ) * 100 };

       return true;
    }
};

class bSphere : public Object
{
    public: double rad2;

            bool type {0};

            Vec3d o_pos{};

            bSphere* owner;

            std::vector<Object*> owns;
            std::vector<Object2*> owns2;

    bSphere( Vec3d p = Vec3d{0,0,0}, double r = 69, bool iv = false, bSphere* own = nullptr, oType ot = oType::none )
    {
        pos = p; rad = r; rad2 = r * r; col = { 42, 36, 24 }; opq = 0.09; isvis = iv;
        owner = own; otype = ot;
        if ( owner ) { owner->bSpheres.push_back(this); }
    }

    virtual void rotate( Vec3d ps, double cx, double sy, int tp ) override

     {   double tx, tz;

         if ( tp == 1 ) { pos = pos - ps; tx=(pos.x*cx)-(pos.z*sy); tz=(pos.z*cx)+(pos.x*sy); pos.z=tz; pos.x=tx; pos = pos + ps; }

         if ( tp == 0 ) { pos = pos - ps; tx=(pos.x*cx)-(pos.y*sy); tz=(pos.y*cx)+(pos.x*sy); pos.y=tz; pos.x=tx; pos = pos + ps; }

         if ( tp == 2 ) { pos = pos - ps; tx=(pos.y*cx)-(pos.z*sy); tz=(pos.z*cx)+(pos.y*sy); pos.z=tz; pos.y=tx; pos = pos + ps; }

         for ( auto o : owns ) { o->rotate(ps,cx,sy,tp); }

         for ( auto b : bSpheres ) { b->rotate(ps,cx,sy,tp); }
    }

   virtual void scale (double scaler) override
   {
     rad2 *= scaler * scaler;
     if ( o_pos.x || o_pos.y || o_pos.z) { pos = ( ( pos - o_pos ) * scaler ) + o_pos; }
   }

   inline virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
   {
       Vec3d L = pos - orig;
       double tca = L.dot(dir);
       if (tca < 0) return false;
       double d2 = L.dot(L) - tca * tca;
       if (d2 > rad2) return false;
       double thc = sqrt( rad2 - d2  );
       dist1 = tca - thc; dist2 = tca + thc;

       dist1 = ( dist1 < dist2 ) ? dist2 : dist1;

       res_col = col;
       return true;
    }
};

class sSphere : public Object
{
    public: double rad2;
            bool type {0};

            Vec3d o_pos{};

            bSphere* owner;


    sSphere( Vec3d p = Vec3d(0,0,0), double r = 69, Vec3d c = Vec3d(0,0,0), double o = 1.0, bSphere* own = nullptr, bool bl = true, bool iv = true )
    {
        pos = p; rad = r; rad2 = r * r; col = c; opq = o; isvis = iv; blur = bl;
        owner = own; otype = oType::sphere;
        if ( owner ) { owner->owns.push_back(this); }
    }

    virtual void rotate(Vec3d, double, double, int ) override
    { }

   virtual void scale (double scaler) override
   {
     rad2 *= scaler * scaler;
     if ( o_pos.x || o_pos.y || o_pos.z) { pos = ( ( pos - o_pos ) * scaler ) + o_pos; }
   }

   inline virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
   {
       Vec3d L = pos - orig;
       double tca = L.dot(dir);
       if (tca < 0) return false;
       double d2 = L.dot(L) - tca * tca;
       if (d2 > rad2) return false;
       double thc = sqrt( rad2 - d2  );
       dist1 = tca - thc; dist2 = tca + thc;

       dist1 = ( dist1 > dist2 ) ? dist2 : dist1;
       Vec3d Phit = ( ( orig + dir * dist1 ) - pos ) . normalize() * 1.1;

       res_col = { col.x + ( Phit.x / M_PI ) * 100.0,
                   col.y - ( Phit.y / M_PI ) * 100.0,
                   col.z + ( Phit.z / M_PI ) * 100.0 };

       return true;
    }
};

class oSphere : public Object
{
    public: double rad2;
            bool type {0};

            Vec3d o_pos{};

            bSphere* owner;

    oSphere( Vec3d p = Vec3d(0.0,0.0,0.0), double r = 69.0, Vec3d c = Vec3d(0.0,0.0,0.0),
             double o = 1.0, bSphere* own = nullptr, bool bl = true, bool iv = true )
    {
        pos = p; rad = r; rad2 = r * r; col = c; opq = o; isvis = iv; blur = bl;
        owner = own; otype = oType::sphere;
        if ( owner ) { owner->owns.push_back(this); }
    }

    virtual void rotate( Vec3d ps, double cx, double sy, int tp ) override

     {
         double tx, tz;

         if ( tp == 1 ) { pos = pos - ps; tx=(pos.x*cx)-(pos.z*sy); tz=(pos.z*cx)+(pos.x*sy); pos.z=tz; pos.x=tx; pos = pos + ps; }

         if ( tp == 0 ) { pos = pos - ps; tx=(pos.x*cx)-(pos.y*sy); tz=(pos.y*cx)+(pos.x*sy); pos.y=tz; pos.x=tx; pos = pos + ps; }

         if ( tp == 2 ) { pos = pos - ps; tx=(pos.y*cx)-(pos.z*sy); tz=(pos.z*cx)+(pos.y*sy); pos.z=tz; pos.y=tx; pos = pos + ps; }
    }


   virtual void scale (double scaler) override
   {
     rad2 *= scaler * scaler;
     if ( o_pos.x || o_pos.y || o_pos.z) { pos = ( ( pos - o_pos ) * scaler ) + o_pos; }
   }

   inline virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
   {

       double rd =0;

       Vec3d L = pos - orig;
       double tca = L.dot(dir);
       if ( tca < 0 ) return false;
       double d2 = L.dot(L) - tca * tca;
       if ( d2 > rad2 - rd ) return false;
       double thc = sqrt( rad2 - rd - d2  );
       dist1 = tca - thc; dist2 = tca + thc;

       dist1 = ( dist1 > dist2 ) ? dist2 : dist1;
       Vec3d Phit = ( ( orig + dir * dist1 ) - pos ) . normalize() * 1.1;

       res_col = { col.x + ( Phit.x / M_PI ) * 100.0,
                   col.y - ( Phit.y / M_PI ) * 100.0,
                   col.z + ( Phit.z / M_PI ) * 100.0  };

       return true;
    }
};

class aTriangle : public Object
{
public: Vec3d a,b,c,n; bSphere* owner; Vec3d ocol, rcol, bor_col; double border, invbord; bool shd{false}, trc{false};

    aTriangle ( Vec3d a0, Vec3d b0, Vec3d c0, double o = 1.0, Vec3d clr = Vec3d(255.0,255.0,255.0),
                bSphere *own = nullptr, bool iv = true, double br = 0.0, Vec3d brcol = Vec3d(13.0,13.0,13.0) )
        {
          a = a0; b = b0; c = c0; opq = o; isvis = iv; pos = a; col = clr; ocol = col;
          owner = own; if (own) { owner->owns.push_back(this); } otype = oType::triangle;
          border = br; invbord = 1.0 - br; bor_col = brcol; n = (a-b).cross(c-b); n.normalize();
        }

    aTriangle ( ) = default;

    virtual  void rotate(Vec3d ps, double cx, double sy, int tp) override

     {   double tx, tz;

         if ( tp == 1 ) { a = a - ps; tx=(a.x*cx)-(a.z*sy); tz=(a.z*cx)+(a.x*sy); a.z=tz; a.x=tx; a = a + ps;
                          b = b - ps; tx=(b.x*cx)-(b.z*sy); tz=(b.z*cx)+(b.x*sy); b.z=tz; b.x=tx; b = b + ps;
                          c = c - ps; tx=(c.x*cx)-(c.z*sy); tz=(c.z*cx)+(c.x*sy); c.z=tz; c.x=tx; c = c + ps;
                          tx=(n.x*cx)-(n.z*sy); tz=(n.z*cx)+(n.x*sy); n.z=tz; n.x=tx; }

         if ( tp == 0 ) { a = a - ps; tx=(a.x*cx)-(a.y*sy); tz=(a.y*cx)+(a.x*sy); a.y=tz; a.x=tx; a = a + ps;
                          b = b - ps; tx=(b.x*cx)-(b.y*sy); tz=(b.y*cx)+(b.x*sy); b.y=tz; b.x=tx; b = b + ps;
                          c = c - ps; tx=(c.x*cx)-(c.y*sy); tz=(c.y*cx)+(c.x*sy); c.y=tz; c.x=tx; c = c + ps;
                          tx=(n.x*cx)-(n.y*sy); tz=(n.y*cx)+(n.x*sy); n.y=tz; n.x=tx; }

         if ( tp == 2 ) { a = a - ps; tx=(a.y*cx)-(a.z*sy); tz=(a.z*cx)+(a.y*sy); a.z=tz; a.y=tx; a = a + ps;
                          b = b - ps; tx=(b.y*cx)-(b.z*sy); tz=(b.z*cx)+(b.y*sy); b.z=tz; b.y=tx; b = b + ps;
                          c = c - ps; tx=(c.y*cx)-(c.z*sy); tz=(c.z*cx)+(c.y*sy); c.z=tz; c.y=tx; c = c + ps;
                          tx=(n.y*cx)-(n.z*sy); tz=(n.z*cx)+(n.y*sy); n.z=tz; n.y=tx; }

    }

    virtual void scale (double scaler) override {         a = ( ( a - owner->pos ) * scaler ) + owner->pos;
                                                          b = ( ( b - owner->pos ) * scaler ) + owner->pos;
                                                          c = ( ( c - owner->pos ) * scaler ) + owner->pos; }

    inline virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
    {
        double u, v, sh;


        Vec3d v0v1 = b - a, v0v2 = c - a,
              pvec = dir.cross(v0v2);

        double det = v0v1.dot(pvec);
        if (det > -kEpsilon && det < kEpsilon) return false;
        double invDet = 1.0 / det;

        Vec3d tvec = orig - a;
        u = tvec.dot(pvec) * invDet;
        if ( u < 0.0 || u > 1.0 ) return false;

        Vec3d qvec = tvec.cross(v0v1);
        v = dir.dot(qvec) * invDet;
        if ( v < 0.0 || v + u > 1.0 ) return false;

        //if ( t < kEpsilon) return false;
        dist1 = v0v2.dot(qvec) * invDet; dist2 = INFINITY;


        sh = { ( std::abs(n.dot(dir))) };

        res_col = Vec3d{ col.x * sh + exp( atan2 (det, u + v ) ), col.y * sh * exp( atan2( cos(v), sin(u) ) / M_E ), col.z * sh * u , };
        return true;

    }

};

class Plane : public Object
{
    public: Vec3d whl0{};

            Vec3d aa{}, ab{}, ac{}, ad{};

            Vec3d p0{};

            Vec3d a1{},  b1{},  c1{},  d1{};

           double  a1l{}, b1l{}, c1l{}, d1l{};

            Vec3d n0{};

            bSphere* owner;

            uchar *dtZ;

            std::vector<Object*> owns;

    Plane ( Vec3d p00 = Vec3d(0.0,0.0,0.0), Vec3d whl = Vec3d(0.0,0.0,0.0), Vec3d clr = Vec3d(35.0,35.0,35.0),
            double o = 1.0, bSphere *own = nullptr, uchar *dtZ0 = nullptr, bool iv = false )
    { pos = p00; whl0 = whl; otype = oType::plane;
      opq = o; col = clr; isvis = iv; dtZ = dtZ0;
      owner = own; if (own) { owner->owns.push_back(this); }

      n0 = { 0.0, 0.0, 1.0};
      p0 = { pos.x, pos.y, pos.z };

          aa = { p0.x + whl0.x / 2.0, p0.y + whl0.y / 2.0, p0.z };
          ab = { p0.x - whl0.x / 2.0, p0.y + whl0.y / 2.0, p0.z };
          ac = { p0.x - whl0.x / 2.0, p0.y - whl0.y / 2.0, p0.z };
          ad = { p0.x + whl0.x / 2.0, p0.y - whl0.y / 2.0, p0.z };

          a1 = aa - ab; b1 = aa - ad; c1 = ac - ad; d1 = ac - ab;

          a1l = a1.length(); b1l = b1.length(); c1l = c1.length(); d1l = d1.length();

    };

    virtual void rotate(Vec3d, double, double, int) override { }

    virtual void scale (double) override { }

    virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
    {
        Vec3d d{dir.x, dir.y, dir.z};

        double cosF1{}, cosF2{}, cosF3{}, cosF4{}, denom = { n0.dot(d) };
        Vec3d ipntv{};

        if ( denom > kEpsilon || denom < -kEpsilon ) {

                             dist1 =  ( p0 - orig ) . dot( n0 ) / denom;
                             Vec3d pnt = ( orig + dir * dist1 );

                             ipntv = pnt - ab;
                             cosF1 = ipntv.dot( a1 ) / ( b1l * ipntv.length() );

                             if ( cosF1 > 0.0 ) {
                             ipntv = pnt - ad;
                             cosF2 = ipntv.dot( c1 ) / ( d1l * ipntv.length() );

                             if ( cosF2 > 0.0 ) {
                             ipntv = pnt - aa;
                             cosF3 = ipntv.dot( d1 ) / ( c1l * ipntv.length() );

                             if ( cosF3 > 0.0 ) {
                             ipntv = pnt - ac;
                             cosF4 = ipntv.dot( b1 ) / ( a1l * ipntv.length() );

                             if ( cosF4 > 0.0 ) { dist1 = ipntv.x; dist2 = ipntv.y;

                             int x{ int( dist1 * 2.0 ) }, y{ int( dist2 * 2.0 ) }, i{ ( ( y * 494 ) + x ) * 3 };

                   res_col = Vec3d {

                             static_cast<double>( dtZ[ i + 0 ] ) * ( cosF3 + cosF4 * sin ( cosF2 + cosF1 ) ) * 0.65 ,
                             static_cast<double>( dtZ[ i + 1 ] ) * ( 1 - sin ( cosF1 + cosF2 ) * 0.75  ) ,
                             static_cast<double>( dtZ[ i + 2 ] ) * ( cos ( cosF4 + cosF3 ) * 0.75 ) };

                             return true; } } } }

           }

        return false;

    }

};

class Line : public Object
{
    public: Vec3d p0{}, p1{}, n0{};

            bSphere* owner;

            std::vector<Object*> owns;

    Line ( Vec3d p00 = Vec3d(0.0,0.0,0.0), Vec3d p01 = Vec3d(0.0,0.0,0.0),
           Vec3d clr = Vec3d(35.0,35.0,35.0), double o = 1.0, bSphere *own = nullptr, bool iv = false )
    { pos = p00; otype = oType::line;
      opq = o; col = clr; isvis = iv;
      owner = own; if (own) { owner->owns.push_back(this); }

      p0 = p00; p1 = p01;
      n0 = ( p1 - p0 ) . normalize();


    };

    virtual void rotate(Vec3d, double, double, int) override { }

    virtual void scale (double) override { }

    virtual bool intersect (const Vec3d &orig, const Vec3d &dir, double &dist1, double &dist2, Vec3d &res_col) override
    {
            Vec3d u = n0 . cross(dir);

            double is = ( orig - p0 ) . dot ( u );

            dist1 = u . norm();

            if ( std::abs( is * is ) < dist1 ) {

                Vec3d pp0 { ( p0 - orig ) . cross(u) }, pp1 { ( p1 - orig ) . cross(u) };

                double esc{ 1.0 }, dd0{ dir . dot(pp0) }, dd1{ dir . dot(pp1) };

                if ( dd0 < 0.0 || dd1 > 0.0 ) return false;

                esc = ( dd0 > 1.0 ) ? ( dd1 < -1.0 ) ? 1.0 : -1.0 * dd1 : dd0;

            //  Vec3d tmp0 = p0 - orig, tmp1 = p1 - orig;
            //  double dp0 { tmp0.dot(dir) }, dp1 { tmp1.dot(dir) };
            //  if ( dp0 <= 0 && dp1 <= 0 ) return false;

                dist1 = std::sqrt(dist1);

                dist2 = ( dist1 - std::abs( is ) ) * esc * opq;

                dist1 = ( ( p0 + ( n0 * dd0 ) ) - orig ) . norm();// std::cout << dist1 << "\n";

                res_col = col;

                return true; }

        return false;

    }

};

class Point : public Object
{
    public: double drad {}, dk{}, deg{};

            bool yn{false};

            bSphere* owner;

            std::vector<Object*> owns;

    Point ( Vec3d p00 = Vec3d(0.0,0.0,0.0), Vec3d clr = Vec3d(35.0,35.0,35.0), double o = 1.0,
            bSphere *own = nullptr, double d0 = 1.0, bool iv = false )
    {
      pos = p00; otype = oType::point;
      opq = o; col = clr; isvis = iv; drad = d0; dk = 1.0 / drad ; opq *= dk;
      owner = own; if (own) { owner->owns.push_back(this); }
    }

    virtual void rotate( Vec3d px, double cx, double sy, int tp ) override
    {
        pos.rotate( px, cx, sy, tp );
    }

    virtual void scale ( double ) override { }

    virtual bool intersect ( const Vec3d &o, const Vec3d &d, double &dist1, double &dist2, Vec3d &res_col ) override
    {
        Vec3d fvec { pos - o },
              pvec { Vec3d ( fvec ) . normalize() },
              tvec { pvec . cross ( d ) };

        dist1 = fvec . norm();

        double pvar { tvec . norm() * dist1 };

        if ( pvar < drad ) { //std::cout << dist1 <<"\n";

                dist2 = opq * ( drad - pvar );

                res_col = col;

                return true;
        }

        return false;

    }

};

class Ellipse : public Object2
{
    public: Vec3d p0{}, n0{}, f0{} ;

            bSphere* owner;

            double scw{8.0}, size{};

            int tp{0};

            std::vector<Object*> owns;

    Ellipse ( Vec3d p00 = Vec3d(0.0,0.0,0.0), bSphere *own = nullptr, double o = 1.0,
              Vec3d clr = Vec3d(35.0,35.0,35.0), double s0 = 500.0, int tp0 = 0, bool iv = false )
    { pos = p00; otype = oType::order2; size = s0; tp = tp0;
      opq = o; col = clr; isvis = iv;
      owner = own; if (own) { owner->owns2.push_back(this); }

      p0 = p00;
      f0 = p0;
      f0.y += 10.0;
      n0.z = 1.0; }

    bool solveQ(double a, double b, double c, double &x0, double &x1)
    {
        double d = b * b - 4.0 * a * c;
        if ( d < 0 ) { return false; }
        else if ( d == 0 ) { x0 = x1 = -0.5 * b / a; }
        else { double q = ( b > 0 ) ? -0.5 * ( b + sqrt( d ) ) : -0.5 * ( b - sqrt( d ) );
               x0 = q / a ; x1 = c / q ; }
        return true;
    }

    virtual void rotate(Vec3d, double, double, int) override { }

    virtual void scale (double) override { }

    virtual bool intersect (const Vec3d &o, const Vec3d &d, double &dist1, double &op1, Vec3d &res_col, double &dist2, double &op2, Vec3d &res_col2) override
    {
        Vec3d fv = o - pos,
              dv{ d.x, d.y / 2, d.z / 3 },
              fv2{ fv.x, fv.y / 2, fv.z / 3 };

        double a = dv.dot(dv),
               b = 2.0 * fv2.dot(dv),
               c = fv2.dot(fv2) - size;

        if ( solveQ( a, b, c, dist1, dist2 ) ) { op1 = op2 = opq;

            //double d1{dist1}, d2{dist2};
            if ( dist1 <= 0 ) { dist1 = INFINITY; } if ( dist2 <= 0 ) { dist2 = INFINITY; }

            if ( dist1 > dist2 ) { std::swap(dist1, dist2); if ( dist1 == INFINITY ) return false; }

            Vec3d Fhit = ( ( o + d * dist1 ) - pos ) . normalize(),
                  Phit = ( ( o + d * dist2 ) - pos ) . normalize();

            double dotd1{ ( dist1 < 0 ) ? d.dot(Fhit) : 1 },
                   dotd2{ ( dist2 < 0 ) ? d.dot(Phit) : 1 };

            uint ais1 {}, ais2 {};

            switch ( tp ) { case 0 : break;

            case 1 : ais1 = uint ( ( ( Fhit.z - atan ( Fhit.x / ( std::abs( Fhit.y ) + 1 ) ) ) ) * scw );
                     ais2 = uint ( ( ( Phit.z - atan ( Phit.x / ( std::abs( Phit.y ) + 1 ) ) ) ) * scw );
                     break;

            case 2 : ais1 = uint ( ( ( Fhit.z - acos ( Fhit.x / atan( Fhit.y ) ) ) ) * scw * 0.5 );
                     ais2 = uint ( ( ( Phit.z - acos ( Phit.x / atan( Phit.y ) ) ) ) * scw * 0.5 );
                     break;

            case 3 : ais1 = uint ( ( ( Fhit.x + atan ( Fhit.z * Fhit.y ) ) ) * scw * 0.5 );
                     ais2 = uint ( ( ( Phit.x + atan ( Phit.z * Phit.y ) ) ) * scw * 0.5 );
                     break;

            case 4 : ais1 = uint ( ( ( Fhit.z * sin ( tan ( Fhit.x ) + asin ( Fhit.y ) ) ) ) * scw );
                     ais2 = uint ( ( ( Phit.z * sin ( tan ( Phit.x ) + asin ( Phit.y ) ) ) ) * scw );
                     break;

            }

            if ( ( tp != 0 && ais1 % 2 == 0 ) || dotd1 < 0 ) { dist1 = INFINITY; } else {

                // double oo{ (double)ais * 0.1 }; if ( oo < 1 && oo > 0) op1 *= oo;

            res_col = { col.x - ( Fhit.x / M_PI ) * 190.0,
                        col.y + ( Fhit.y / M_PI ) * 190.0,
                        col.z + ( Fhit.z / M_PI ) * 190.0 }; }

            if ( ( tp != 0 && ais2 % 2 == 0 ) || dotd2 < 0 ) { dist2 = INFINITY; } else {

                // double oo{ (double)ais * 0.1 }; if ( oo < 1 && oo > 0) op2 *= oo;

            res_col2 = { col.x - ( Phit.x / M_PI ) * 190.0,
                         col.y + ( Phit.y / M_PI ) * 190.0,
                         col.z + ( Phit.z / M_PI ) * 190.0 }; }

       return true; } return false;

    }

};

class Hyper1 : public Object2
    {
        public: Vec3d p0{}, n0{};

                double sci{1.0}, scj{1.0}, sck{1.0},
                       s1i{1.0}, s1j{1.0}, s1k{1.0},
                       s2i{1.0}, s2j{1.0}, s2k{1.0},

                       scc{100.0}, scw{8.0}, angx{0.0}, size{500.0};

                bool sect{false};
                     bool sepr{false};

                Matrix44d rmx;

                bSphere *owner;

                char *dtZ{nullptr};

                std::vector<Object*> owns;

        Hyper1 ( Vec3d p00 = Vec3d(0.0,0.0,0.0), bSphere *own = nullptr, Vec3d clr = Vec3d(45.0,55.0,45.0),
                 bool iv = false )
        { pos = p00; otype = oType::order2;
          col = clr; isvis = iv;
          owner = own; if (own) { owner->owns2.push_back(this); }

          p0 = p00; opq = 1.0;

        };

        bool solveQ(double a, double b, double c, double &x0, double &x1)
        {
            double d = b * b - 4.0 * a * c;
            if ( d < 0 ) { return false; }
            else if ( d == 0 ) { x0 = x1 = -0.5 * b / a; }
            else { double q = ( b > 0 ) ? -0.5 * ( b + sqrt( d ) ) : -0.5 * ( b - sqrt( d ) );
                   x0 = q / a ; x1 = c / q ; }
            return true;
        }

        virtual void rotate(Vec3d, double, double, int) override { }

        virtual void scale (double) override { }

        virtual bool intersect (const Vec3d &o, const Vec3d &d, double &dist1, double &op1, Vec3d &res_col, double &dist2, double &op2, Vec3d &res_col2) override
        {
            Vec3d fv = o - pos,
                  tdv  = { d.x / s1i, d.y / s1j, d.z / s1k },
                  tfv2 = { fv.x / s1i , fv.y / s1j, fv.z / s1k },
                  dv  = {}, dv2 = {},
                  fv2 = {}, fv3 = {};

            rmx.multDirMatrix( tdv, dv );
            rmx.multDirMatrix( tfv2, fv2 );

            dv2 = { sci * dv.x / s2i, scj * dv.y / s2j, sck * dv.z / s2k };
            fv3 = { sci * fv2.x / s2i, scj * fv2.y / s2j, sck * fv2.z / s2k };

            double a = dv.dot(dv2);
            double b = 2.0 * fv2.dot(dv2);
            double c = fv3.dot(fv2) - size;

            if ( solveQ( a, b, c, dist1, dist2 ) ) { double sh{1};

                double d1{dist1}, d2{dist2};
                if ( dist1 <= 0 ) { dist1 = INFINITY; } if ( dist2 <= 0 ) { dist2 = INFINITY; }

                if ( dist1 > dist2 ) { dist1 = dist2; std::swap(d1, d2); }
                if ( dist1 == INFINITY ) return false;

                Vec3d Phit = ( ( o + d * dist1 ) - pos ), Thit{Phit}; rmx.multDirMatrix( Phit, Thit );

                double opq1 = ( Phit.norm() - 22500.0 ) / 2500.0;
                       opq1 = ( opq1 < 0 ) ? 1.0 : ( opq1 > 1.0 ) ? 0.0 : 1.0 - opq1;

                int xx{ (int)Thit.x * 2 + 680 },yy{ (int)Thit.y * 2 + 680 };

                Vec3d addc{};

                if ( dtZ && xx >= 0 && yy >= 0 && xx < 1360 && yy < 1360 ) { int ps = ( ( yy * 1360 ) + xx ) * 3;
                addc = ( Vec3d{ (double)(static_cast<uchar>(dtZ[ps+0])), (double)(static_cast<uchar>(dtZ[ps+1])), (double)(static_cast<uchar>(dtZ[ps+2])) } ) * 0.16; sh = 1;
                } else { addc = Vec3d{0,0,0}; sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize() ) * 0.5 ); }

                Phit . normalize(); if ( sect ) {
                                                  if ( std::abs(Phit.x) < 0.3 ) { dist1 = INFINITY; goto ex1_; }
                                                  if ( std::abs(Phit.x) > 0.3 && ( std::abs(Phit.x) < 0.4 ) )
                                                     { double o0 = ( std::abs(Phit.x) - 0.3 ) * 10; opq1 *= o0; }
                                                }
                op1 = opq1 * opq;
                                    if ( sepr ) {
                                                  uint ais1 =

                                                  uint ( ( ( sin ( cos ( Thit.z ) * sin ( Thit.y ) ) + sin ( cos ( Thit.y ) + sin ( Thit.x ) ) ) ) * scw );

                                                  if ( ais1 % 4 == 0 ) { dist1 = INFINITY; goto ex1_; }
                                                }



                   res_col = Vec3d { col.x + ( Phit.x / M_PI ) * 100.0 ,
                                     col.y - ( Phit.y / M_PI ) * 100.0 ,
                                     col.z + ( Phit.z / M_PI ) * 100.0 } * sh + addc ;
ex1_:
                    double opq2; dist2 = d2;

                    Phit = ( ( o + d * ( d2 ) ) - pos ); rmx.multDirMatrix( Phit, Thit );

                           opq2 = ( Phit.norm() - 22500.0 ) / 2500.0;
                           opq2 = ( opq2 < 0.0 ) ? 1.0 : ( opq2 > 1.0 ) ? 0.0 : 1.0 - opq2;

                    xx = (int)Thit.x * 2 + 680.0; yy = (int)Thit.y * 2 + 680.0;

                    if ( dtZ && xx >= 0 && yy >= 0  && xx < 1360 && yy < 1360 ) { int ps = ( ( yy * 1360 ) + xx ) * 3;
                    addc = ( Vec3d{ (double)(static_cast<uchar>(dtZ[ps+0])), (double)(static_cast<uchar>(dtZ[ps+1])), (double)(static_cast<uchar>(dtZ[ps+2])) } ) * 0.16; sh = 1;
                    } else { addc = Vec3d{0,0,0}; sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize() ) * 0.5 ); }

                    Phit . normalize(); if ( sect ) {
                                                      if ( std::abs(Phit.x) < 0.3 ) dist2 = INFINITY;
                                                      if ( std::abs(Phit.x) > 0.3 && ( std::abs(Phit.x) < 0.4 ) )
                                                         { double o0 = ( std::abs(Phit.x) - 0.3 ) * 10; opq2 *= o0; }
                                                   }
                                        if ( sepr ) {
                                                            uint ais1 =

                                                            uint ( ( ( sin ( cos ( Thit.z ) * sin ( Thit.y ) ) + sin ( cos ( Thit.y ) + sin ( Thit.x ) ) ) ) * scw );

                                                            if ( ais1 % 4 == 0 ) { dist2 = INFINITY;
                                                                if ( dist1 == dist2 ) return false; else return true; }
                                                         }

                    op2 = opq2 * opq;

                    res_col2 = Vec3d { col.x + ( Phit.x / M_PI ) * 100.0 ,
                                       col.y - ( Phit.y / M_PI ) * 100.0 ,
                                       col.z + ( Phit.z / M_PI ) * 100.0 } * sh + addc;

           return true; } return false;

        }

};

class Hyper2 : public Object2
    {
        public: Vec3d p0{}, n0{};

                bSphere* owner;

                double sci{1.0}, scj{1.0}, sck{1.0},
                       s1i{1.0}, s1j{1.0}, s1k{1.0},
                       s2i{1.0}, s2j{1.0}, s2k{1.0},

                       scc{100.0}, angx{90};

                double scw{8.0};

                bool sect{false}, sepr{false};

                char *dtZ{nullptr};

                Matrix44d rmx;

                std::vector<Object*> owns;

        Hyper2 ( Vec3d p00 = Vec3d(0.0,0.0,0.0), bSphere *own = nullptr, Vec3d clr = Vec3d(45.0,55.0,45.0),
                 bool iv = false )
        { pos = p00; otype = oType::order2;
          col = clr; isvis = iv;
          owner = own; if (own) { owner->owns2.push_back(this); }
          p0 = p00; };

        inline bool solveQ(double a, double b, double c, double &x0, double &x1)
        {
            double d = b * b - 4.0 * a * c;
            if ( d < 0 ) { return false; }
            else if ( d == 0 ) { x0 = x1 = -0.5 * b / a; }
            else { double q = ( b > 0 ) ? -0.5 * ( b + sqrt( d ) ) : -0.5 * ( b - sqrt( d ) );
                   x0 = q / a ; x1 = c / q ; }
            return true;
        }

        virtual void rotate(Vec3d, double, double, int) override { }

        virtual void scale (double) override { }

        virtual bool intersect (const Vec3d &o, const Vec3d &d, double &dist1, double &op1, Vec3d &res_col, double &dist2, double &op2, Vec3d &res_col2) override
        {
            Vec3d fv = o - pos,
                  tdv = { d.x / s1i, d.y / s1j, d.z / s1k },
                  tfv2 = { fv.x / s1i , fv.y / s1j, fv.z / s1k },
                  dv  = {}, dv2 = {},
                  fv2 = {}, fv3 = {};

            rmx.multDirMatrix(tdv,dv);
            rmx.multDirMatrix(tfv2,fv2);

            dv2 ={ sci * dv.x / s2i, scj * dv.y / s2j, sck * dv.z / s2k };
            fv3 ={ sci * fv2.x / s2i, scj * fv2.y / s2j, sck * fv2.z / s2k };

            double a = dv.dot(dv2);
            double b = 2.0 * fv2.dot(dv2);
            //if ( b > 0 ) return false;
            double c = fv3.dot(fv2) - 250;

            if ( solveQ( a, b, c, dist1, dist2 ) ) { double sh{1.0};

                double d1{dist1}, d2{dist2};
                if ( dist1 <= 0 ) { dist1 = INFINITY; } if ( dist2 <= 0 ) { dist2 = INFINITY; }

                if ( dist1 > dist2 ) { dist1 = dist2; std::swap(d1, d2); }
                if ( dist1 == INFINITY ) return false;

                Vec3d Phit = ( ( o + d * dist1 ) - pos ), Thit{}; rmx.multDirMatrix( Phit, Thit );
                Vec3d Zhit = Thit;

                double opq = ( Phit.norm() - 2500.0 ) / 22500.0;
                       opq = ( opq < 0 ) ? 1.0 : ( opq > 1.0 ) ? 0.0 : 1.0 - opq;

                       int xx{ (int)Thit.z * 2 + 680 },yy{ (int)Thit.y * 2 + 680 };

                       Vec3d addc{};

                       if ( dtZ && xx >= 0 && yy >= 0 && xx < 1360 && yy < 1360 ) { int ps = ( ( yy * 1360 ) + xx ) * 3;
                       addc = ( Vec3d{ (double)(static_cast<uchar>(dtZ[ps+0])), (double)(static_cast<uchar>(dtZ[ps+1])), (double)(static_cast<uchar>(dtZ[ps+2])) } ) * 0.16; sh = 1;
                       } else { addc = Vec3d{0,0,0}; sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize() ) * 0.5 ); }

                Phit . normalize(); if ( sect ) { if ( std::abs(Phit.z) < 0.5 ) dist1 = INFINITY;
                                                  if ( std::abs(Phit.z) > 0.5 && ( std::abs(Phit.z) < 0.6 ) )
                                                     { double o0 = ( std::abs(Phit.z) - 0.5 ) * 10; opq *= o0; }
                                                }
                                    else if ( sepr ) {
                                                        uint ais1 =
                                                        uint ( ( ( Zhit.x ) ) ); //- ( cos ( Thit.x ) + sin ( Thit.y ) ) ) ) * scw ); //* scw
                                                        if ( ais1 % 3 == 0 ) goto ex1_;
                                                     }

                op1 = opq;

                sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize()) * 0.5 );

                res_col = Vec3d { col.x + ( Phit.x / M_PI ) * scc ,
                                  col.y - ( Phit.y / M_PI ) * scc ,
                                  col.z + ( Phit.z / M_PI ) * scc } * sh + addc;

        ex1_:
                    double opq2; dist2 = d2;

                    Phit = ( ( o + d * ( d2 ) ) - pos ); rmx.multDirMatrix( Phit, Thit ); Zhit = Thit;

                           opq2 = ( Phit.norm() - 2500.0 ) / 22500.0;
                           opq2 = ( opq2 < 0.0 ) ? 1.0 : ( opq2 > 1.0 ) ? 0.0 : 1.0 - opq2;

                    xx = (int)Thit.z * 2 + 680.0; yy = (int)Thit.y * 2 + 680.0;

                    if ( dtZ && xx >= 0 && yy >= 0  && xx < 1360 && yy < 1360 ) { int ps = ( ( yy * 1360 ) + xx ) * 3;
                    addc = ( Vec3d{ (double)(static_cast<uchar>(dtZ[ps+0])), (double)(static_cast<uchar>(dtZ[ps+1])), (double)(static_cast<uchar>(dtZ[ps+2])) } ) * 0.16; sh = 1;
                    } else { addc = Vec3d{0,0,0}; sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize() ) * 0.5 ); }

                    Phit . normalize(); if ( sect ) { if ( std::abs(Phit.z) < 0.5 ) dist2 = INFINITY;
                                                      if ( std::abs(Phit.z) > 0.5 && ( std::abs(Phit.z) < 0.6 ) )
                                                         { double o0 = ( std::abs(Phit.z) - 0.5 ) * 10; opq2 *= o0; }
                                                    }
                                        else if ( sepr ) {
                                                            uint ais1 =
                                                            uint ( ( ( Zhit.x ) ) ); // - ( cos ( Thit.x ) + sin ( Thit.y ) ) ) ) * scw ); //* scw
                                                            if ( ais1 % 3 == 0 ) return false;
                                                         }

                    op2 = opq2;

                    sh = 1 + std::abs( 0 - d.dot( ( Phit - pos ).normalize()) * 0.5 );

                    res_col2 = Vec3d { col.x + ( Phit.x / M_PI ) * scc ,
                                       col.y - ( Phit.y / M_PI ) * scc ,
                                       col.z + ( Phit.z / M_PI ) * scc } * sh;

           return true; } return false;

        }

};




#endif // GEOMETRY_H
