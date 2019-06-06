//
// Created by Henrique on 05/06/2019.
//

#ifndef RAYTRACERCHALLENGE_MATRIX_H
#define RAYTRACERCHALLENGE_MATRIX_H

#include "tuples.h"
#include "util.h"

struct Matrix2x2{
    public:
        float m[2][2] = { 0 };

        Matrix2x2(){
            m[0][0] = 1; m[0][1] = 0;
            m[1][0] = 0; m[1][1] = 1;
        }

        Matrix2x2(float m00, float m01,
                  float m10, float m11){

            m[0][0] = m00; m[0][1] = m01;
            m[1][0] = m10; m[1][1] = m11;

        }

        float* operator [](unsigned int pos);
};

struct Matrix3x3{
    public:
        float m[3][3] = { 0 };

        Matrix3x3(){
            m[0][0] = 1; m[0][1] = 0; m[0][2] = 0;
            m[1][0] = 0; m[1][1] = 1; m[1][2] = 0;
            m[2][0] = 0; m[2][1] = 0; m[2][2] = 1;
        }

        Matrix3x3(float m00, float m01, float m02,
                  float m10, float m11, float m12,
                  float m20, float m21, float m22){

            m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
            m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
            m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;

        }

        float* operator [](unsigned int pos);
};

struct Matrix4x4{
    public:
        float m[4][4] = { 0 };

        Matrix4x4(){

            m[0][0] = 1; m[0][1] = 0; m[0][2] = 0; m[0][3] = 0;
            m[1][0] = 0; m[1][1] = 1; m[1][2] = 0; m[1][3] = 0;
            m[2][0] = 0; m[2][1] = 0; m[2][2] = 1; m[2][3] = 0;
            m[3][0] = 0; m[3][1] = 0; m[3][2] = 0; m[3][3] = 1;

        }

        Matrix4x4(float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33){

            m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
            m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
            m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
            m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;

        }

        float* operator [](unsigned int pos);
};

float* Matrix2x2::operator [](unsigned int pos){
    return m[pos];
}

float* Matrix3x3::operator [](unsigned int pos){
    return m[pos];
}

float* Matrix4x4::operator [](unsigned int pos){
    return m[pos];
}

bool operator == (Matrix4x4 m1, Matrix4x4 m2){
    for(unsigned int i = 0; i < 4; i++){
        for(unsigned int j = 0; j < 4; j++){
            if(!equal(m1[i][j], m2[i][j])){
                return false;
            }
        }
    }
    return true;
}

bool operator != (Matrix4x4 m1, Matrix4x4 m2){
    return !(m1 == m2);
}

bool operator == (Matrix3x3 m1, Matrix3x3 m2){
    for(unsigned int i = 0; i < 3; i++){
        for(unsigned int j = 0; j < 3; j++){
            if(!equal(m1[i][j], m2[i][j])){
                return false;
            }
        }
    }
    return true;
}

bool operator != (Matrix3x3 m1, Matrix3x3 m2){
    return !(m1 == m2);
}

bool operator == (Matrix2x2 m1, Matrix2x2 m2){
    for(unsigned int i = 0; i < 2; i++){
        for(unsigned int j = 0; j < 2; j++){
            if(!equal(m1[i][j], m2[i][j])){
                return false;
            }
        }
    }
    return true;
}

bool operator != (Matrix2x2 m1, Matrix2x2 m2){
    return !(m1 == m2);
}

Matrix4x4 multiply(Matrix4x4 m1, Matrix4x4 m2){
    auto resultado = Matrix4x4();

    for(int row = 0; row < 4; row++){
        for(int col = 0; col < 4; col++){
            resultado[row][col] = (m1[row][0] * m2[0][col]) +
                                  (m1[row][1] * m2[1][col]) +
                                  (m1[row][2] * m2[2][col]) +
                                  (m1[row][3] * m2[3][col]);
        }
    }

    return resultado;
}

Matrix4x4 operator * (Matrix4x4 m1, Matrix4x4 m2){
    return multiply(m1, m2);
}

Tuple multiply(Matrix4x4 m, Tuple t){
    auto resultado = Tuple{};

    resultado.x = (t.x * m[0][0]) + (t.y * m[0][1]) + (t.z * m[0][2]) + (t.w * m[0][3]);
    resultado.y = (t.x * m[1][0]) + (t.y * m[1][1]) + (t.z * m[1][2]) + (t.w * m[1][3]);
    resultado.z = (t.x * m[2][0]) + (t.y * m[2][1]) + (t.z * m[2][2]) + (t.w * m[2][3]);
    resultado.w = (t.x * m[3][0]) + (t.y * m[3][1]) + (t.z * m[3][2]) + (t.w * m[3][3]);

    return resultado;
}

Tuple operator * (Matrix4x4 m, Tuple t){
    return multiply(m, t);
}

Tuple operator * (Tuple t, Matrix4x4 m){
    return multiply(m, t);
}

Matrix4x4 identity(){
    return {};
}

Matrix4x4 transpose(Matrix4x4 m){
    auto transpose = Matrix4x4();

    for(int x = 0; x < 4; x++){
        for(int y = 0; y < 4; y++){
            transpose[x][y] = m[y][x];
        }
    }

    return transpose;
}

float determinant(Matrix2x2 m){
    return (m[0][0] * m[1][1]) - (m[0][1] * m[1][0]);
}

Matrix3x3 submatrix(Matrix4x4 m, unsigned int row, unsigned int column){

    auto sub = Matrix3x3();

    if(row < 4 && column < 4){
        int rowSub = 0, colSub = 0;
        for(int x = 0; x < 4; x++){
            if(x != row){
                for(int y = 0; y < 4; y++){
                    if(y != column){
                        sub[rowSub][colSub++] = m[x][y];
                    }
                }
                rowSub++;
                colSub = 0;
            }
        }
    }

    return sub;
}

Matrix2x2 submatrix(Matrix3x3 m, unsigned int row, unsigned int column){

    auto sub = Matrix2x2();

    if(row < 3 && column < 3){
        int rowSub = 0, colSub = 0;
        for(int x = 0; x < 3; x++){
            if(x != row){
                for(int y = 0; y < 3; y++){
                    if(y != column){
                        sub[rowSub][colSub++] = m[x][y];
                    }
                }
                rowSub++;
                colSub = 0;
            }
        }
    }

    return sub;
}

float minor(Matrix3x3 m, unsigned int row, unsigned int column){
    return determinant(submatrix(m, row, column));
}

float cofactor(Matrix3x3 m, unsigned int row, unsigned int column){
    auto cf = minor(m, row, column);
    if((row + column) % 2 != 0){
        cf = -cf;
    }
    return cf;
}

float determinant(Matrix3x3 m){
    auto cf1 = cofactor(m, 0, 0);
    auto cf2 = cofactor(m, 0, 1);
    auto cf3 = cofactor(m, 0, 2);
    auto det = (m[0][0] * cf1) + (m[0][1] * cf2) + (m[0][2] * cf3);
    return det;
}

float minor(Matrix4x4 m, unsigned int row, unsigned int column){
    return determinant(submatrix(m, row, column));
}

float cofactor(Matrix4x4 m, unsigned int row, unsigned int column){
    auto cf = minor(m, row, column);
    if((row + column) % 2 != 0){
        cf = -cf;
    }
    return cf;
}

float determinant(Matrix4x4 m){
    return (m[0][0] * cofactor(m, 0, 0)) + (m[0][1] * cofactor(m, 0, 1)) + (m[0][2] * cofactor(m, 0, 2)) + (m[0][3] * cofactor(m, 0, 3));
}

Matrix4x4 inverse(Matrix4x4 m){
    auto inv = Matrix4x4();
    auto det = determinant(m);

    //Check if m is invertible
    if(det != 0){
        for(unsigned int row = 0; row < 4; row++){
            for(unsigned int col = 0; col < 4; col++){
                auto cf = cofactor(m, row, col);
                inv[col][row] = cf / det;
            }
        }
    }

    return inv;
}

void printMatrix(Matrix4x4 m){
    for(int x = 0; x < 4; x++){
        for(int y = 0; y < 4; y++){
            printf("%f", m[x][y]);
            if(y < 3){
                printf(", ");
            }
        }
        printf("\n");
    }
}

void printMatrix(Matrix3x3 m){
    for(int x = 0; x < 3; x++){
        for(int y = 0; y < 3; y++){
            printf("%f", m[x][y]);
            if(y < 2){
                printf(", ");
            }
        }
        printf("\n");
    }
}

void printMatrix(Matrix2x2 m){
    for(int x = 0; x < 2; x++){
        for(int y = 0; y < 2; y++){
            printf("%f", m[x][y]);
            if(y < 1){
                printf(", ");
            }
        }
        printf("\n");
    }
}

#endif //RAYTRACERCHALLENGE_MATRIX_H
