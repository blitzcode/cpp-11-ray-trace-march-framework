
#include "lin_alg.h"

void TestLinAlg()
{
    // Instantiate templates, check sizes, make sure operators compile etc.

    static_assert(sizeof(Vec2f)  == 8 , "Vec2f size test failed" );
    static_assert(sizeof(Vec3f)  == 12, "Vec3f size test failed" );
    static_assert(sizeof(Vec4f)  == 16, "Vec4f size test failed" );
    static_assert(sizeof(Vec2d)  == 16, "Vec2d size test failed" );
    static_assert(sizeof(Vec3d)  == 24, "Vec3d size test failed" );
    static_assert(sizeof(Vec4d)  == 32, "Vec4d size test failed" );
    static_assert(sizeof(Vec2i)  == 8 , "Vec2i size test failed" );
    static_assert(sizeof(Vec3i)  == 12, "Vec3i size test failed" );
    static_assert(sizeof(Vec4i)  == 16, "Vec4i size test failed" );
    static_assert(sizeof(Vec2ui) == 8 , "Vec2ui size test failed");
    static_assert(sizeof(Vec3ui) == 12, "Vec3ui size test failed");
    static_assert(sizeof(Vec4ui) == 16, "Vec4ui size test failed");

    Vec2f vec2f(0.0f, 0.0f);
    Vec3f vec3f(0.0f, 0.0f, 0.0f);
    Vec4f vec4f(0.0f, 0.0f, 0.0f, 0.0f);
    Vec2d vec2d(0.0, 0.0);
    Vec3d vec3d(0.0, 0.0, 0.0);
    Vec4d vec4d(0.0, 0.0, 0.0, 0.0);
    Vec2i vec2i(-1, -1);
    Vec3i vec3i(-1, -1, -1);
    Vec4i vec4i(-1, -1, -1, -1);
    Vec2ui vec2ui(0, 0);
    Vec3ui vec3ui(0, 0, 0);
    Vec4ui vec4ui(0, 0, 0, 0);
    float f = 0.0f;
    bool b = false;

    b = (vec3f == vec3f);
    b = (vec3f != vec3f);
    vec3f = Vec3f(1.0f) + Vec3f(2.0f);
    vec3f = Vec3f(1.0f) - Vec3f(2.0f);
    vec3f = Vec3f(1.0f) * Vec3f(2.0f);
    vec3f = Vec3f(1.0f) / Vec3f(2.0f);
    vec3f = Vec3f(1.0f) * f;
    vec3f = f * Vec3f(1.0f);
    vec3f = Vec3f(1.0f) / f;
    vec3f = -Vec3f(1.0f);
    vec3f += Vec3f(1.0f);
    vec3f -= Vec3f(1.0f);
    vec3f *= Vec3f(1.0f);
    vec3f /= Vec3f(1.0f);
    vec3f *= f;
    vec3f /= f;
    f = vec3f[0];
    vec3f[0] = f;
    f = Length(vec3f);
    f = LengthSquared(vec3f);
    vec3f = Normalize(vec3f);
    f = Dot(Vec3f(1.0f), Vec3f(2.0f));
    vec3f = Vec3f(1.0f) ^ Vec3f(2.0f);

    Matrix44f matf;
    matf.Identity();
    Matrix44d matd;
    matf.RotationX(1);
    matf.RotationY(1);
    matf.RotationZ(1);
    b = matf == matf;
    matf = matf * matf;
    matf.BuildLookAtMatrix(Vec3f(0.0f, 10.0f, 10.0f), Vec3f(0.0f));
    matf.BuildProjection(90.0f, 4.0f / 3.0f, 1.0f, 1000.0f);
    Vec3f out;
    matf.Transf4x4(vec3f, out);
    matf.Transpose3x3();
    matf.Transpose4x4();
    matf.Invert();
}

