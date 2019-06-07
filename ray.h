//
// Created by Henrique on 06/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RAY_H
#define RAYTRACERCHALLENGE_RAY_H

#include <cstdlib>
#include "tuples.h"
#include "matrix.h"
#include "canvas.h"


struct Light{
    Tuple position;
    Color intensity;
};

struct Ray {
    Tuple origin;
    Tuple direction;
};

Ray ray(float originX, float originY, float originZ, float directionX, float directionY, float directionZ){
    return Ray {point(originX, originY, originZ), vector(directionX, directionY, directionZ)};
}

Tuple position(Ray ray, float time){
    return ray.origin + (ray.direction * time);
}

Ray transform(Ray ray, Matrix4x4 matrix){
    auto newRay = Ray{};
    newRay.origin = ray.origin * matrix;
    newRay.direction = ray.direction * matrix;
    return newRay;
}

Tuple reflect(Tuple v, Tuple normal){
    return v - (normal * 2 * dot(v, normal));
}

Light pointLight(Tuple position, Color intensity){
    return Light{point(position.x, position.y, position.z), intensity};
}

#endif //RAYTRACERCHALLENGE_RAY_H
