#include <iostream>
#include <cmath>
#include "Helpers.h"

using namespace std;

void setVec3(Vec3 & a, const Vec3 & b)
{
    a.x = b.x;
    a.y = b.y;
    a.z = b.z;
}

Vec3 multiplyVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.y * b.z - b.y * a.z, b.x * a.z - a.x * b.z, a.x * b.y - b.x * a.y);
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v)
{
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v)
{
    double d = magnitudeOfVec3(v);
    return Vec3(v.x / d, v.y / d, v.z / d);
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v)
{
    return Vec3(-v.x, -v.y, -v.z);
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b)
{
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c)
{
    return Vec3(v.x * c, v.y * c, v.z * c);
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v)
{
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << std::endl;
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b)
{

    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((ABS((a.x - b.x)) < EPSILON) && (ABS((a.y - b.y)) < EPSILON) && (ABS((a.z - b.z)) < EPSILON))
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
 */
Matrix4 getIdentityMatrix()
{
    Matrix4 result;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (i == j)
            {
                result.values[i][j] = 1.0;
            }
            else
            {
                result.values[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.values[i][k] * m2.values[k][j];
            }

            result.values[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = 0;
        for (int j = 0; j < 4; j++)
        {
            total += m.values[i][j] * v.getNthComponent(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

Vec3 find_v(const Vec3 &u) {
    if (u.x < u.y) 
    {
        if (u.x < u.z) return Vec3(0, u.z, -u.y);
        else return Vec3(u.y, -u.x, 0); //u.z < u.x
    }
    else //u.y < u.x
    {
        if (u.y < u.z) return Vec3(u.z, 0, -u.x);
        else return Vec3(u.y, -u.x, 0);//u.z < u.y
    }
}

Vec3 rotateVec3(const Vec3 &vec, const Vec3 &u, double angle) {

    Vec3 v = find_v(u);
    v = normalizeVec3(v);

    Vec3 w = crossProductVec3(u, v);
    w = normalizeVec3(w);

    Matrix4 M = Matrix4({
        {u.x, u.y, u.z},
        {v.x, v.y, v.z},
        {w.x, w.y, w.z},
        {0,0,0,1}
    });

    Matrix4 M_inv = Matrix4({
        {u.x, v.x, w.x, 0},
        {u.y, v.y, w.y, 0},
        {u.z, v.z, w.z, 0},
        {0,0,0,1}
    });

    Matrix4 Rx_angle = Matrix4({
        {1,0,0,0},
        {0, cos(angle), -sin(angle), 0},
        {0, sin(angle), cos(angle), 0},
        {0,0,0,1}
    });

    Matrix4 rot_matrix = multiplyMatrixWithMatrix(M_inv, multiplyMatrixWithMatrix(Rx_angle, M)); 

    Vec4 vec4 = Vec4(vec.x, vec.y, vec.z, 1, vec.colorId);

    vec4 = multiplyMatrixWithVec4(rot_matrix, vec4);

    return Vec3(vec4.x, vec4.y, vec4.z, vec4.colorId);

    //rodrigues
    // return addVec3(
    //     addVec3(multiplyVec3WithScalar(vec, cos(angle)), 
    //             multiplyVec3WithScalar(crossProductVec3(u, vec), sin(angle))),
    //     multiplyVec3WithScalar(u, dotProductVec3(u, vec) * (1 - cos(angle)))
    // );

}

void swap(Color & c0, Color & c1) {
    swap(c0.r, c1.r);
    swap(c0.g, c1.g);
    swap(c0.b, c1.b);
}

int min3(int a, int b, int c) {
    if (a < b) 
    {
        if (a < c) return a;
    }
    else 
    {
        if (b < c) return b;
    }
    return c;
}

int max3(int a, int b, int c) {
    if (a > b) 
    {
        if (a > c) return a;
    }
    else 
    {
        if (b > c) return b;
    }
    return c;
}

double line_eqn(const Vec4 & vec0, const Vec4 & vec1, double x, double y) {
    return x*(vec0.y - vec1.y) + y*(vec1.x - vec0.x) + vec0.x*vec1.y - vec0.y*vec1.x;
}