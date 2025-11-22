#ifndef __MATRIX4_H__
#define __MATRIX4_H__
#include <iostream>

class Matrix4
{
public:
    double values[4][4];

    Matrix4();
    Matrix4(std::initializer_list<std::initializer_list<double>> init) {
        int i = 0, j;
        for (std::initializer_list<double> row : init) {
            j = 0;
            for (double val : row) {
                values[i][j] = val;
                j++;
            }
            ++i;
        }
    }
    Matrix4(double values[4][4]);
    Matrix4(const Matrix4 &other);
    friend std::ostream &operator<<(std::ostream &os, const Matrix4 &m);
};

#endif