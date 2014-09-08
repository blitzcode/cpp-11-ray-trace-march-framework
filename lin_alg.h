
#ifndef LIN_ALG
#define LIN_ALG

#include <algorithm>
#include <cmath>
#include <cstring>
#include <type_traits>

#include "types.h"

//
// Vector
//

// Storage, explicit specializations for 2, 3 and 4 components
template<typename T, uint N> struct VectorStorage_t
{
    T m_vec[N]; // General case, no named components
};
template<typename T> struct VectorStorage_t<T, 2>
{
    union { struct { T x; T y; }; T m_vec[2]; };
};
template<typename T> struct VectorStorage_t<T, 3>
{
    union { struct { T x; T y; T z; }; T m_vec[3]; };
};
template<typename T> struct VectorStorage_t<T, 4>
{
    union { struct { T x; T y; T z; T w; }; T m_vec[4]; };
};

// Basic operations that every vector type has
template<typename T, uint N> struct BasicVector_t : public VectorStorage_t<T, N>
{
    // Array operator
    const T& operator [] (uint i) const  { return this->m_vec[i]; }
    T& operator       [] (uint i)        { return this->m_vec[i]; }

    // Equality
    bool operator == (const VectorStorage_t<T, N>& a) const
    {
        for (uint i=0; i<N; i++)
            if (this->m_vec[i] != a.m_vec[i])
                return false;
        return true;
    }
    bool operator != (const VectorStorage_t<T, N>& a) const
    {
        for (uint i=0; i<N; i++)
            if (this->m_vec[i] != a.m_vec[i])
                return true;
        return false;
    }
};

// Main vector template. You don't get much if you don't hit one of the specializations
template<typename T, uint N> struct Vector_t : public BasicVector_t <T, N> { };

// Operator definitions which require Vector_t
template<typename T, uint N> struct Operator_t : public BasicVector_t <T, N>
{
    #define IMPL_OP(op_)             \
        Vector_t<T, N> result;       \
        for (uint i=0; i<N; i++)     \
            result.m_vec[i] = (op_); \
        return result;

    // Global operators
    friend Vector_t<T, N> operator * (T scalar, const Vector_t<T, N>& vec)
        { IMPL_OP(vec.m_vec[i] * scalar); }

    // Arithmetic operators
    Vector_t<T, N> operator + (const Vector_t<T, N> &a) const
        { IMPL_OP(this->m_vec[i] + a.m_vec[i]); }
    Vector_t<T, N> operator - (const Vector_t<T, N> &a) const
        { IMPL_OP(this->m_vec[i] - a.m_vec[i]); }
    Vector_t<T, N> operator * (const Vector_t<T, N> &a) const
        { IMPL_OP(this->m_vec[i] * a.m_vec[i]); }
    Vector_t<T, N> operator / (const Vector_t<T, N> &a) const
        { IMPL_OP(this->m_vec[i] / a.m_vec[i]); }
    Vector_t<T, N> operator + (T a) const { IMPL_OP(this->m_vec[i] + a); }
    Vector_t<T, N> operator - (T a) const { IMPL_OP(this->m_vec[i] - a); }
    Vector_t<T, N> operator * (T a) const { IMPL_OP(this->m_vec[i] * a); }
    Vector_t<T, N> operator / (T a) const { IMPL_OP(this->m_vec[i] / a); }
    Vector_t<T, N> operator - ()    const { IMPL_OP(-this->m_vec[i]);    }
    void operator += (const Vector_t<T, N> &a)
        { for (uint i=0; i<N; i++) this->m_vec[i] += a.m_vec[i]; }
    void operator -= (const Vector_t<T, N> &a)
        { for (uint i=0; i<N; i++) this->m_vec[i] -= a.m_vec[i]; }
    void operator *= (const Vector_t<T, N> &a)
        { for (uint i=0; i<N; i++) this->m_vec[i] *= a.m_vec[i]; }
    void operator /= (const Vector_t<T, N> &a)
        { for (uint i=0; i<N; i++) this->m_vec[i] /= a.m_vec[i]; }
    void operator *= (T a)
        { for (uint i=0; i<N; i++) this->m_vec[i] *= a; }
    void operator /= (T a)
        { for (uint i=0; i<N; i++) this->m_vec[i] /= a; }

    #undef IMPL_OP
};

// Operations only defined for three dimensions
template<typename T, uint N> struct Dim3Only_t : public Operator_t<T, N> { };
template<typename T> struct Dim3Only_t<T, 3> : public Operator_t<T, 3>
{
    friend Vector_t<T, 3> operator ^ (const Vector_t<T, 3>& a, const Vector_t<T, 3>& b)
    {
        return Vector_t<T, 3>
            (a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    friend Vector_t<T, 3> Cross(const Vector_t<T, 3>& a, const Vector_t<T, 3>& b)
    {
        return a ^ b;
    }
    friend Vector_t<T, 3> Reflect(const Vector_t<T, 3>& n, const Vector_t<T, 3>& i)
    {
        // R = I - 2 * (N*I) * N
        // Note that this assumes I is pointing towards the plane of reflection
        static_assert(std::is_floating_point<T>(), "needs to be a floating point type");
        T angle = Dot(n, i) * T(2.0);
        return i - n * angle;
    }
    friend uint32 ToBGRA8(const Vector_t<T, 3>& col)
    {
        static_assert(std::is_floating_point<T>(), "needs to be a floating point type");
        uchar rc = col.x > 1.0f ? 255 : (uchar) (col.x * 255.0f);
        uchar gc = col.y > 1.0f ? 255 : (uchar) (col.y * 255.0f);
        uchar bc = col.z > 1.0f ? 255 : (uchar) (col.z * 255.0f);
        return rc << 16 | gc << 8  | bc << 0;
    }
};

// Common base class containing global functions operating on vectors
template<typename T, uint N> struct CommonBase_t : public Dim3Only_t<T, N>
{
    friend T Dot(const Vector_t<T, N>& a, const Vector_t<T, N>& b)
    {
        T result = T();
        for (uint i=0; i<N; i++)
            result += a.m_vec[i] * b.m_vec[i];
        return result;
    }
    friend T LengthSquared(const Vector_t<T, N>& vec) { return Dot(vec, vec); }
    friend T Length(const Vector_t<T, N>& vec) { return std::sqrt(LengthSquared(vec)); }
    friend T DistanceSquared(const Vector_t<T, N>& a, const Vector_t<T, N>& b)
        { return LengthSquared(a - b); }
    friend T Distance(const Vector_t<T, N>& a, const Vector_t<T, N>& b)
        { return std::sqrt(DistanceSquared(a, b)); }
    friend Vector_t<T, N> Normalize(const Vector_t<T, N>& vec)
    {
        static_assert(std::is_floating_point<T>(), "needs to be a floating point type");
        const T length = T(1.0) / Length(vec);
        return vec * length;
    }
    friend Vector_t<T, N> ComponentMin(const Vector_t<T, N>& vec1, const Vector_t<T, N>& vec2)
    {
        Vector_t<T, N> result;
        for (uint i=0; i<N; i++)
            result[i] = std::min(vec1[i], vec2[i]);
        return result;
    }
    friend Vector_t<T, N> ComponentMax(const Vector_t<T, N>& vec1, const Vector_t<T, N>& vec2)
    {
        Vector_t<T, N> result;
        for (uint i=0; i<N; i++)
            result[i] = std::max(vec1[i], vec2[i]);
        return result;
    }
};

// Vector class with construction specializations for 2, 3 and 4 components
template<typename T> struct Vector_t<T, 2> : public CommonBase_t<T, 2>
{
    Vector_t()                       { }
    Vector_t(T x_, T y_)             { Set(x_, y_); }
    explicit Vector_t(T scalar)      { this->x = scalar; this->y = scalar; }
    explicit Vector_t(const T *vec)  { this->x = vec[0]; this->y = vec[1]; }
    void Set(T x_, T y_)             { this->x = x_; this->y = y_; }
};
template<typename T> struct Vector_t<T, 3> : public CommonBase_t<T, 3>
{
    Vector_t()                       { }
    Vector_t(T x_, T y_, T z_)       { Set(x_, y_, z_); }
    explicit Vector_t(T scalar)      { this->x = scalar; this->y = scalar; this->z = scalar; }
    explicit Vector_t(const T *vec)  { this->x = vec[0]; this->y = vec[1]; this->z = vec[2]; }
    void Set(T x_, T y_, T z_)       { this->x = x_; this->y = y_; this->z = z_; }
};
template<typename T> struct Vector_t<T, 4> : public CommonBase_t<T, 4>
{
    Vector_t() { }
    Vector_t(T x_, T y_, T z_, T w_) { Set(x_, y_, z_, w_); }
    explicit Vector_t(T scalar)      { this->x = scalar;
                                       this->y = scalar;
                                       this->z = scalar;
                                       this->w = scalar; }
    explicit Vector_t(const T *vec)  { this->x = vec[0];
                                       this->y = vec[1];
                                       this->z = vec[2];
                                       this->w = vec[3]; }
    void Set(T x_, T y_, T z_, T w_) { this->x = x_; this->y = y_; this->z = z_; this->w = w_; }
};

template <class T> T Clamp(T val, T min, T max)
{
    if (val < min)
        return min;
    if (val > max)
        return max;
    return val;
}

// Common type / size definitions
typedef Vector_t<float,  2> Vec2f;
typedef Vector_t<float,  3> Vec3f;
typedef Vector_t<float,  4> Vec4f;
typedef Vector_t<double, 2> Vec2d;
typedef Vector_t<double, 3> Vec3d;
typedef Vector_t<double, 4> Vec4d;
typedef Vector_t<int,    2> Vec2i;
typedef Vector_t<int,    3> Vec3i;
typedef Vector_t<int,    4> Vec4i;
typedef Vector_t<uint,   2> Vec2ui;
typedef Vector_t<uint,   3> Vec3ui;
typedef Vector_t<uint,   4> Vec4ui;

//
// Matrix
//

template<typename T> T DegToRad(const T deg) { return deg * T(0.0174532925 ); }
template<typename T> T RadToDeg(const T rad) { return rad * T(57.2957795131); }

template <class T> struct Matrix44_t
{
    Matrix44_t()            { Identity();                            }
    Matrix44_t(T mat[4][4]) { std::memcpy(mat, mat, sizeof(T) * 16); }

    Matrix44_t(T f11, T f12, T f13,
               T f21, T f22, T f23,
               T f31, T f32, T f33)
    {
        Set(f11, f12, f13,
            f21, f22, f23,
            f31, f32, f33);
    }

    Matrix44_t(T f11, T f12, T f13, T f14,
               T f21, T f22, T f23, T f24,
               T f31, T f32, T f33, T f34,
               T f41, T f42, T f43, T f44)
    {
        Set(f11, f12, f13, f14,
            f21, f22, f23, f24,
            f31, f32, f33, f34,
            f41, f42, f43, f44);
    }

    void Set(T f11, T f12, T f13,
             T f21, T f22, T f23,
             T f31, T f32, T f33)
    {
        Identity();

        m_mat[0][0] = f11; m_mat[1][0] = f12; m_mat[2][0] = f13;
        m_mat[0][1] = f21; m_mat[1][1] = f22; m_mat[2][1] = f23;
        m_mat[0][2] = f31; m_mat[1][2] = f32; m_mat[2][2] = f33;
    }

    void Set(T f11, T f12, T f13, T f14,
             T f21, T f22, T f23, T f24,
             T f31, T f32, T f33, T f34,
             T f41, T f42, T f43, T f44)
    {
        m_mat[0][0] = f11; m_mat[1][0] = f12; m_mat[2][0] = f13; m_mat[3][0] = f14;
        m_mat[0][1] = f21; m_mat[1][1] = f22; m_mat[2][1] = f23; m_mat[3][1] = f24;
        m_mat[0][2] = f31; m_mat[1][2] = f32; m_mat[2][2] = f33; m_mat[3][2] = f34;
        m_mat[0][3] = f41; m_mat[1][3] = f42; m_mat[2][3] = f43; m_mat[3][3] = f44;
    }

    bool operator == (const Matrix44_t& a) const
    {
        uint i, j;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
            {
                if (m_mat[i][j] != a.m_mat[i][j])
                    return false;
            }
        return true;
    }

    bool operator != (const Matrix44_t& a) const
    {
        uint i, j;
        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
            {
                if (m_mat[i][j] != a.m_mat[i][j])
                    return true;
            }
        return false;
    }

    friend Matrix44_t operator * (const Matrix44_t& a, const Matrix44_t& b)
        { Matrix44_t mat_ret(a); mat_ret.Multiply(b); return mat_ret; };

    void CopyTo(T mat[4][4]) const { std::memcpy(mat, m_mat, sizeof(T) * 16); };

    void Identity()
    {
        std::memset(m_mat, 0, sizeof(T) * 16);
        m_mat[0][0] = m_mat[1][1] = m_mat[2][2] = m_mat[3][3] = T(1.0);
    }

    void Translation(T x, T y, T z)
    {
        // |  1    0    0    tx |   x' = x + tx
        // |  0    1    0    ty |   y' = y + ty
        // |  0    0    1    tz |   z' = z + tz
        // |  0    0    0    1  |

        Set(T(1.0), T(0.0), T(0.0), x,
            T(0.0), T(1.0), T(0.0), y,
            T(0.0), T(0.0), T(1.0), z,
            T(0.0), T(0.0), T(0.0), T(1.0));
    }

    void RotationX(T degrees)
    {
        T r = DegToRad(degrees);

        // |  1      0      0     0 |   x' = x
        // |  0    cos a  -sin a  0 |   y' = (cos a) * y - (sin a) * z
        // |  0    sin a   cos a  0 |   z' = (sin a) * y + (cos a) * z
        // |  0      0      0     1 |

        Set(T(1.0), T(0.0),         T(0.0),          T(0.0),
            T(0.0), T(std::cos(r)), T(-std::sin(r)), T(0.0),
            T(0.0), T(std::sin(r)), T(std::cos(r)),  T(0.0),
            T(0.0), T(0.0),         T(0.0),          T(1.0));
    }

    void RotationY(T degrees)
    {
        T r = DegToRad(degrees);

        // |  cos a   0    -sin a   0 |  x' = (cos a) * x - (sin a) * z
        // |    0     1     0       0 |  y' = y
        // |  sin a   0    cos a    0 |  z' = (sin a) * x + (cos a) * z
        // |    0     0     0       1 |

        Set(T(std::cos(r)), T(0.0), T(-std::sin(r)), T(0.0),
            T(0.0),         T(1.0), T(0.0),          T(0.0),
            T(std::sin(r)), T(0.0), T(std::cos(r)),  T(0.0),
            T(0.0),         T(0.0), T(0.0),          T(1.0));
    }

    void RotationZ(T degrees)
    {
        T r = DegToRad(degrees);

        // |  cos a   -sin a   0    0 | x' = (cos a) * x - (sin a) * y
        // |  sin a    cos a   0    0 | y' = (sin a) * x + (cos a) * y
        // |    0       0      1    0 | z' = z
        // |    0       0      0    1 |

        Set(T(std::cos(r)), T(-std::sin(r)), T(0.0), T(0.0),
            T(std::sin(r)), T(std::cos(r)),  T(0.0), T(0.0),
            T(0.0),         T(0.0),          T(1.0), T(0.0),
            T(0.0),         T(0.0),          T(0.0), T(1.0));
    }

    void BuildBasis3x3(const Vector_t<T, 3>& normal,
                       const Vector_t<T, 3>& binormal,
                       const Vector_t<T, 3>& tangent)
    {
        Set(normal.x,   normal.y,   normal.z,
            binormal.x, binormal.y, binormal.z,
            tangent.x,  tangent.y,  tangent.z);
    }

    void AxisRotation(const Vector_t<T, 3>& axis, T degrees)
    {
        // Compute a matrix which rotates degrees around axis. We do this by
        // finding a basis for axis, transforming into this basis, rotating and
        // transforming back. See Realtime Rendering 2nd Edition, page 42

        //      t
        // X = M * R (a) * M
        //          x

        Vector_t<T, 3> binormal, tangent;
        Matrix44_t mat_transform, mat_transform_t, mat_rot;

        // Create M and MT
        axis.ComputeBasis(binormal, tangent);
        mat_transform.BuildBasis3x3(axis, binormal, tangent);
        mat_transform_t = mat_transform;
        mat_transform_t.Transpose3x3();

        // Rotate around new x
        mat_rot.RotationX(degrees);

        // Build final transform
        (* this) = mat_transform_t * mat_rot * mat_transform;
    }

    void BuildProjection(T fov, T aspect, T near, T far)
    {
        T f43, f11, f22;
        T left, right, top, bottom;

        bottom  = -(near * T(std::tan(DegToRad(fov / T(2.0)))));
        left    = bottom;
        bottom /= aspect;
        top     = -bottom;
        right   = -left;

        f11 = -1.0f * near / (right - left);
        f22 = -1.0f * near / (top - bottom);
        f43 = 2.0f * far * near / (far - near);

        Set(f11,    T(0.0), T(0.0),  T(0.0),
            T(0.0), f22,    T(0.0),  T(0.0),
            T(0.0), T(0.0), T(1.0),  f43,
            T(0.0), T(0.0), T(-1.0), T(0.0));
    }

    void BuildLookAtMatrix(const Vector_t<T, 3>& eye,
                           const Vector_t<T, 3>& look_at,
                           const Vector_t<T, 3>& up = Vector_t<T, 3>(T(0.0), T(1.0), T(0.0)))
    {
        Matrix44_t mat_look_at_transf;
        mat_look_at_transf.Identity();

        // Z Axis points at the object
        Vector_t<T, 3> zaxis = Normalize(eye - look_at);

        // X Axis perpendicular to Up and Z
        Vector_t<T, 3> xaxis = Normalize(up ^ zaxis);

        // Y Axis perpendicular to Z and X
        Vector_t<T, 3> yaxis = Normalize(zaxis ^ xaxis);

        // Build basis
        mat_look_at_transf.m_mat[0][0] = xaxis.x;
        mat_look_at_transf.m_mat[0][1] = xaxis.y;
        mat_look_at_transf.m_mat[0][2] = xaxis.z;
        mat_look_at_transf.m_mat[1][0] = yaxis.x;
        mat_look_at_transf.m_mat[1][1] = yaxis.y;
        mat_look_at_transf.m_mat[1][2] = yaxis.z;
        mat_look_at_transf.m_mat[2][0] = zaxis.x;
        mat_look_at_transf.m_mat[2][1] = zaxis.y;
        mat_look_at_transf.m_mat[2][2] = zaxis.z;

        // Compute translation
        Matrix44_t mat_trans;
        mat_trans.Translation(
            Dot(xaxis, eye),
            Dot(yaxis, eye),
            Dot(zaxis, eye));

        // Combine translation + rotation
        * this = mat_trans * mat_look_at_transf;
    }

    void Scaling(T factor)
    {
        set(factor, T(0.0), T(0.0), T(0.0),
            T(0.0), factor, T(0.0), T(0.0),
            T(0.0), T(0.0), factor, T(0.0),
            T(0.0), T(0.0), T(0.0), T(1.0));
    }

    void Multiply(const Matrix44_t& a)
    {
        uint i, j;
        T mat_old[4][4];

        std::memcpy(mat_old, m_mat, sizeof(T) * 16);

        for (i=0; i<4; i++)
            for (j=0; j<4; j++)
            {
                m_mat[i][j] =
                    mat_old[i][0] * a.m_mat[0][j] +
                    mat_old[i][1] * a.m_mat[1][j] +
                    mat_old[i][2] * a.m_mat[2][j] +
                    mat_old[i][3] * a.m_mat[3][j];
            }
    };

    void Transf3x3(const Vector_t<T, 3>& point, Vector_t<T, 3>& point_out) const
    {
        point_out.x =
            point.x * m_mat[0][0] +
            point.y * m_mat[1][0] +
            point.z * m_mat[2][0];
        point_out.y =
            point.x * m_mat[0][1] +
            point.y * m_mat[1][1] +
            point.z * m_mat[2][1];
        point_out.z =
            point.x * m_mat[0][2] +
            point.y * m_mat[1][2] +
            point.z * m_mat[2][2];
    }

    void Transf3x3(Vector_t<T, 3>& point) const
    {
        Vector_t<T, 3> t;
        Transf3x3(point, t);
        point = t;
    }

    void Transf4x4(const Vector_t<T, 3>& point, Vector_t<T, 3>& point_out) const
    {
        point_out.x =
            point.x * m_mat[0][0] +
            point.y * m_mat[1][0] +
            point.z * m_mat[2][0] +
            m_mat[3][0];
        point_out.y =
            point.x * m_mat[0][1] +
            point.y * m_mat[1][1] +
            point.z * m_mat[2][1] +
            m_mat[3][1];
        point_out.z =
            point.x * m_mat[0][2] +
            point.y * m_mat[1][2] +
            point.z * m_mat[2][2] +
            m_mat[3][2];
    }

    void Transf4x4(Vector_t<T, 3>& point) const
    {
        Vector_t<T, 3> t;
        Transf4x4(point, t);
        point = t;
    }

    void Transf4x4Homogenous(const Vector_t<T, 3> &v, Vector_t<T, 3> &out) const
    {
        T w;

        // Normal 3x3 transform
        out.x = v.x * m_mat[0][0] +
                v.y * m_mat[1][0] +
                v.z * m_mat[2][0] +
                      m_mat[3][0];

        out.y = v.x * m_mat[0][1] +
                v.y * m_mat[1][1] +
                v.z * m_mat[2][1] +
                      m_mat[3][1];

        out.z = v.x * m_mat[0][2] +
                v.y * m_mat[1][2] +
                v.z * m_mat[2][2] +
                      m_mat[3][2];

        // Homogeneous coordinate
        w    = v.x * m_mat[0][3] +
               v.y * m_mat[1][3] +
               v.z * m_mat[2][3] +
                     m_mat[3][3];

        // Avoid division by zero
        assert(w != T(0.0));

        // Convert from homogeneous coordinates to parallel
        out /= w;
    }

    void Transf4x4Homogenous(Vector_t<T, 3> &v) const
    {
        Vector_t<T, 3> t;
        Transf4x4Homogenous(v, t);
        v = t;
    }

    void Transpose4x4()
    {
        Matrix44_t mat_temp = * this;

        m_mat[0][0] = mat_temp.m_mat[0][0];
        m_mat[0][1] = mat_temp.m_mat[1][0];
        m_mat[0][2] = mat_temp.m_mat[2][0];
        m_mat[0][3] = mat_temp.m_mat[3][0];

        m_mat[1][0] = mat_temp.m_mat[0][1];
        m_mat[1][1] = mat_temp.m_mat[1][1];
        m_mat[1][2] = mat_temp.m_mat[2][1];
        m_mat[1][3] = mat_temp.m_mat[3][1];

        m_mat[2][0] = mat_temp.m_mat[0][2];
        m_mat[2][1] = mat_temp.m_mat[1][2];
        m_mat[2][2] = mat_temp.m_mat[2][2];
        m_mat[2][3] = mat_temp.m_mat[3][2];

        m_mat[3][0] = mat_temp.m_mat[0][3];
        m_mat[3][1] = mat_temp.m_mat[1][3];
        m_mat[3][2] = mat_temp.m_mat[2][3];
        m_mat[3][3] = mat_temp.m_mat[3][3];
    }

    void Transpose3x3()
    {
        Matrix44_t mat_temp = * this;

        m_mat[0][0] = mat_temp.m_mat[0][0];
        m_mat[0][1] = mat_temp.m_mat[1][0];
        m_mat[0][2] = mat_temp.m_mat[2][0];

        m_mat[1][0] = mat_temp.m_mat[0][1];
        m_mat[1][1] = mat_temp.m_mat[1][1];
        m_mat[1][2] = mat_temp.m_mat[2][1];

        m_mat[2][0] = mat_temp.m_mat[0][2];
        m_mat[2][1] = mat_temp.m_mat[1][2];
        m_mat[2][2] = mat_temp.m_mat[2][2];
    }

    void Transpose3x3NegTranslation()
    {
        Transpose3x3();

        m_mat[3][0] = -m_mat[3][0];
        m_mat[3][1] = -m_mat[3][1];
        m_mat[3][2] = -m_mat[3][2];
    }

    bool Invert()
    {
        // From MESA's GLU implementation
        // http://cgit.freedesktop.org/mesa/glu/tree/src/libutil/project.c

        const T *m = &m_mat[0][0];
        T inv[16];

        inv[0 ] =  m[5] * m[10] * m[15] - m[5 ] * m[11] * m[14] - m[9 ] * m[6 ] * m[15] +
                   m[9] * m[7 ] * m[14] + m[13] * m[6 ] * m[11] - m[13] * m[7 ] * m[10];
        inv[4 ] = -m[4] * m[10] * m[15] + m[4 ] * m[11] * m[14] + m[8 ] * m[6 ] * m[15] -
                   m[8] * m[7 ] * m[14] - m[12] * m[6 ] * m[11] + m[12] * m[7 ] * m[10];
        inv[8 ] =  m[4] * m[9 ] * m[15] - m[4 ] * m[11] * m[13] - m[8 ] * m[5 ] * m[15] +
                   m[8] * m[7 ] * m[13] + m[12] * m[5 ] * m[11] - m[12] * m[7 ] * m[9 ];
        inv[12] = -m[4] * m[9 ] * m[14] + m[4 ] * m[10] * m[13] + m[8 ] * m[5 ] * m[14] -
                   m[8] * m[6 ] * m[13] - m[12] * m[5 ] * m[10] + m[12] * m[6 ] * m[9 ];
        inv[1 ] = -m[1] * m[10] * m[15] + m[1 ] * m[11] * m[14] + m[9 ] * m[2 ] * m[15] -
                   m[9] * m[3 ] * m[14] - m[13] * m[2 ] * m[11] + m[13] * m[3 ] * m[10];
        inv[5 ] =  m[0] * m[10] * m[15] - m[0 ] * m[11] * m[14] - m[8 ] * m[2 ] * m[15] +
                   m[8] * m[3 ] * m[14] + m[12] * m[2 ] * m[11] - m[12] * m[3 ] * m[10];
        inv[9 ] = -m[0] * m[9 ] * m[15] + m[0 ] * m[11] * m[13] + m[8 ] * m[1 ] * m[15] -
                   m[8] * m[3 ] * m[13] - m[12] * m[1 ] * m[11] + m[12] * m[3 ] * m[9 ];
        inv[13] =  m[0] * m[9 ] * m[14] - m[0 ] * m[10] * m[13] - m[8 ] * m[1 ] * m[14] +
                   m[8] * m[2 ] * m[13] + m[12] * m[1 ] * m[10] - m[12] * m[2 ] * m[9 ];
        inv[2 ] =  m[1] * m[6 ] * m[15] - m[1 ] * m[7 ] * m[14] - m[5 ] * m[2 ] * m[15] +
                   m[5] * m[3 ] * m[14] + m[13] * m[2 ] * m[7 ] - m[13] * m[3 ] * m[6 ];
        inv[6 ] = -m[0] * m[6 ] * m[15] + m[0 ] * m[7 ] * m[14] + m[4 ] * m[2 ] * m[15] -
                   m[4] * m[3 ] * m[14] - m[12] * m[2 ] * m[7 ] + m[12] * m[3 ] * m[6 ];
        inv[10] =  m[0] * m[5 ] * m[15] - m[0 ] * m[7 ] * m[13] - m[4 ] * m[1 ] * m[15] +
                   m[4] * m[3 ] * m[13] + m[12] * m[1 ] * m[7 ] - m[12] * m[3 ] * m[5 ];
        inv[14] = -m[0] * m[5 ] * m[14] + m[0 ] * m[6 ] * m[13] + m[4 ] * m[1 ] * m[14] -
                   m[4] * m[2 ] * m[13] - m[12] * m[1 ] * m[6 ] + m[12] * m[2 ] * m[5 ];
        inv[3 ] = -m[1] * m[6 ] * m[11] + m[1 ] * m[7 ] * m[10] + m[5 ] * m[2 ] * m[11] -
                   m[5] * m[3 ] * m[10] - m[9 ] * m[2 ] * m[7 ] + m[9 ] * m[3 ] * m[6 ];
        inv[7 ] =  m[0] * m[6 ] * m[11] - m[0 ] * m[7 ] * m[10] - m[4 ] * m[2 ] * m[11] +
                   m[4] * m[3 ] * m[10] + m[8 ] * m[2 ] * m[7 ] - m[8 ] * m[3 ] * m[6 ];
        inv[11] = -m[0] * m[5 ] * m[11] + m[0 ] * m[7 ] * m[9 ] + m[4 ] * m[1 ] * m[11] -
                   m[4] * m[3 ] * m[9 ] - m[8 ] * m[1 ] * m[7 ] + m[8 ] * m[3 ] * m[5 ];
        inv[15] =  m[0] * m[5 ] * m[10] - m[0 ] * m[6 ] * m[9 ] - m[4 ] * m[1 ] * m[10] +
                   m[4] * m[2 ] * m[9 ] + m[8 ] * m[1 ] * m[6 ] - m[8 ] * m[2 ] * m[5 ];

        T det =
            m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
        if (det == T(0.0))
            return false;

        det = T(1.0) / det;

        for (uint i=0; i<16; i++)
            (&m_mat[0][0])[i] = inv[i] * det;

        return true;
    }

    T m_mat[4][4];
};

typedef Matrix44_t<float> Matrix44f;
typedef Matrix44_t<double> Matrix44d;

#endif // LIN_ALG

