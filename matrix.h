//
// Created by Henrique on 05/06/2019.
//

#ifndef RAYTRACERCHALLENGE_MATRIX_H
#define RAYTRACERCHALLENGE_MATRIX_H

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

#endif //RAYTRACERCHALLENGE_MATRIX_H
