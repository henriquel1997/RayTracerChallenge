//
// Created by Henrique on 04/06/2019.
//

#ifndef RAYTRACERCHALLENGE_UTIL_H
#define RAYTRACERCHALLENGE_UTIL_H

#include "float.h"
#include <limits>

#define EPSILON 0.00001
#define PI 3.14159265358979323846f
#define MAXDOUBLE DBL_MAX
#define MINDOUBLE -DBL_MAX

double absolute(double f){
    if(f < 0){
        f = f * -1;
    }
    return f;
}

bool equal(double f1, double f2){
    return absolute(f1 - f2) < EPSILON;
}

double radians(double degrees){
    return (degrees / 180) * PI;
}

double clamp(double val, double min, double max){
    if(val > max){
        return max;
    }
    if(val < min){
        return min;
    }
    return val;
}

double max(double v1, double v2){
    if(v1 > v2){
        return v1;
    }
    return v2;
}

double min(double v1, double v2){
    if(v1 < v2){
        return v1;
    }
    return v2;
}

double max(double v1, double v2, double v3){
    return max(max(v1, v2), v3);
}

double min(double v1, double v2, double v3){
    return min(min(v1, v2), v3);
}

#endif //RAYTRACERCHALLENGE_UTIL_H
