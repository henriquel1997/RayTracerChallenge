//
// Created by Henrique on 06/06/2019.
//

#ifndef RAYTRACERCHALLENGE_TRANSFORMATIONS_H
#define RAYTRACERCHALLENGE_TRANSFORMATIONS_H

#include "matrix.h"

Matrix4x4 translation(double x, double y, double z){
    auto id = Matrix4x4();
    id[0][3] = x;
    id[1][3] = y;
    id[2][3] = z;
    return id;
}

Matrix4x4 scaling(double x, double y, double z){
    auto id = Matrix4x4();
    id[0][0] = x;
    id[1][1] = y;
    id[2][2] = z;
    return id;
}

Matrix4x4 rotationX(double angle){
    auto id = Matrix4x4();
    auto seno = sin(angle);
    auto cosseno = cos(angle);
    id[1][1] = cosseno;
    id[1][2] = -seno;
    id[2][1] = seno;
    id[2][2] = cosseno;
    return id;
}

Matrix4x4 rotationY(double angle){
    auto id = Matrix4x4();
    auto seno = sin(angle);
    auto cosseno = cos(angle);
    id[0][0] = cosseno;
    id[0][2] = seno;
    id[2][0] = -seno;
    id[2][2] = cosseno;
    return id;
}

Matrix4x4 rotationZ(double angle){
    auto id = Matrix4x4();
    auto seno = sin(angle);
    auto cosseno = cos(angle);
    id[0][0] = cosseno;
    id[0][1] = -seno;
    id[1][0] = seno;
    id[1][1] = cosseno;
    return id;
}

Matrix4x4 shearing(double xy, double xz, double yx, double yz, double zx, double zy){
    auto id = Matrix4x4();
    id[0][1] = xy;
    id[0][2] = xz;
    id[1][0] = yx;
    id[1][2] = yz;
    id[2][1] = zx;
    id[2][2] = zy;
    return id;
}

#endif //RAYTRACERCHALLENGE_TRANSFORMATIONS_H
