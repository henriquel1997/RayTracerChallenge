//
// Created by Henrique on 06/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RAY_H
#define RAYTRACERCHALLENGE_RAY_H

#include <cstdlib>
#include <vector>
#include "tuples.h"
#include "matrix.h"
#include "canvas.h"

struct Ray {
    Tuple origin;
    Tuple direction;
};

struct Sphere{
    unsigned int id;
    Tuple origin;
    float radius;
    Matrix4x4 transform;
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
    return Sphere{cont++, point(0, 0, 0), 1, Matrix4x4()};
}

Ray transform(Ray ray, Matrix4x4 matrix){
    auto newRay = Ray{};
    newRay.origin = ray.origin * matrix;
    newRay.direction = ray.direction * matrix;
    return newRay;
}

std::vector<Intersection> intersect(Ray ray, Sphere sphere){
    auto tRay = transform(ray, inverse(sphere.transform));

    auto sphere_to_ray = tRay.origin - sphere.origin;

    auto a = dot(tRay.direction, tRay.direction);
    auto b = 2 * dot(tRay.direction, sphere_to_ray);
    auto c = dot(sphere_to_ray, sphere_to_ray) - 1;

    auto discriminant = (b * b) - 4 * a * c;

    if(discriminant < 0){
        return std::vector<Intersection>();
    }

    auto raiz = sqrtf(discriminant);
    auto divisor = (2 * a);

    auto t1 = (-b - raiz) / divisor;
    auto t2 = (-b + raiz) / divisor;

    auto lista = std::vector<Intersection>();

    lista.push_back(Intersection{sphere.id, t1});
    lista.push_back(Intersection{sphere.id, t2});

    return lista;
}

Intersection hit(const std::vector<Intersection>& intersections){
    auto hit = Intersection{0, -1};
    for(auto intersection: intersections){
        if(intersection.time >= 0 && (intersection.time < hit.time || hit.time == -1)){
            hit = intersection;
        }
    }
    return hit;
}

void castRays(){

    unsigned int canvasPixels = 1000;

    auto c = canvas(canvasPixels, canvasPixels);
    auto s = sphere();

    auto rayOrigin = point(0, 0, -5);
    float wallSize = 7;
    float wallZ = 10;
    float pixelSize = wallSize / canvasPixels;
    float half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto position = point(worldX, worldY, wallZ);
            auto ray = Ray{rayOrigin, normalize(position - rayOrigin)};

            if(hit(intersect(ray, s)).time >= 0){
                writePixel(&c, x, y, Color{1, 0, 0});
            }
        }
    }

    canvasToPNG(&c, "sphere.png");
}

#endif //RAYTRACERCHALLENGE_RAY_H
