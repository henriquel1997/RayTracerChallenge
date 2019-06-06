//
// Created by Henrique on 04/06/2019.
//

#ifndef RAYTRACERCHALLENGE_UTIL_H
#define RAYTRACERCHALLENGE_UTIL_H

#define EPSILON 0.00001
#define PI 3.14159265358979323846f

float absolute(float f){
    if(f < 0){
        f = f * -1;
    }
    return f;
}

bool equal(float f1, float f2){
    return absolute(f1 - f2) < EPSILON;
}

float radians(float degrees){
    return (degrees / 180) * PI;
}

#endif //RAYTRACERCHALLENGE_UTIL_H
