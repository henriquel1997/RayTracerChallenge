//
// Created by Henrique on 06/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RAY_H
#define RAYTRACERCHALLENGE_RAY_H

#include <cstdlib>
#include "tuples.h"
#include <vector>

struct Ray {
    Tuple origin;
    Tuple direction;
};

struct Sphere{
    unsigned int id;
    Tuple origin;
    float radius;
};

struct Intersection{
    unsigned int id;
    float time;
};

Ray ray(float originX, float originY, float originZ, float directionX, float directionY, float directionZ){
    return Ray {point(originX, originY, originZ), vector(directionX, directionY, directionZ)};
}

Tuple position(Ray* ray, float time){
    return ray ->origin + (ray->direction * time);
}

Sphere sphere(){
    static unsigned int cont = 0;
    return Sphere{cont++, point(0, 0, 0), 1};
}

std::vector<Intersection> intersect(Ray* ray, Sphere* sphere){
    auto sphere_to_ray = ray->origin - sphere->origin;

    auto a = dot(&ray->direction, &ray->direction);
    auto b = 2 * dot(&ray->direction, &sphere_to_ray);
    auto c = dot(&sphere_to_ray, &sphere_to_ray) - 1;

    auto discriminant = (b * b) - 4 * a * c;

    if(discriminant < 0){
        return std::vector<Intersection>();
    }

    auto raiz = sqrtf(discriminant);
    auto divisor = (2 * a);

    auto t1 = (-b - raiz) / divisor;
    auto t2 = (-b + raiz) / divisor;

    auto lista = std::vector<Intersection>();

    lista.push_back(Intersection{sphere->id, t1});
    lista.push_back(Intersection{sphere->id, t2});

    return lista;
}

#endif //RAYTRACERCHALLENGE_RAY_H
