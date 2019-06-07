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

struct Intersection{
    unsigned int id;
    float time;
};

struct Light{
    Tuple position;
    Color intensity;
};

struct Material{
    Color color;
    float ambient;
    float diffuse;
    float specular;
    float shininess;
};

struct Ray {
    Tuple origin;
    Tuple direction;
};

struct Sphere{
    unsigned int id;
    Tuple origin;
    float radius;
    Matrix4x4 transform;
    Material material;
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

Material material(){
    return Material{ Color{1, 1, 1}, 0.1f, 0.9f, 0.9f, 200 };
}

Sphere sphere(){
    static unsigned int cont = 0;
    return Sphere{cont++, point(0, 0, 0), 1, Matrix4x4(), material()};
}

Tuple normalAt(Sphere sphere, Tuple p){
    auto objectPoint = inverse(sphere.transform) * p;
    auto objectNormal = objectPoint - sphere.origin;
    auto worldNormal = transpose(inverse(sphere.transform)) * objectNormal;
    worldNormal.w = 0;
    return normalize(worldNormal);
}

Tuple reflect(Tuple v, Tuple normal){
    return v - (normal * 2 * dot(v, normal));
}

Light pointLight(Tuple position, Color intensity){
    return Light{point(position.x, position.y, position.z), intensity};
}

Color lighting(Material material, Light light, Tuple point, Tuple eyev, Tuple normalv){
    //Combine the surface color with the light's color/intensity
    auto effectiveColor = material.color * light.intensity;

    //Find the direction to the light source
    auto lightv = normalize(light.position - point);

    //Compute the ambient contribution
    auto ambient = effectiveColor * material.ambient;

    //lightDotNormal represents the cosine of the angle between the light vector and the normal vector.
    //A negative number means the light is on the other side of the surface.
    auto lightDotNormal = dot(lightv, normalv);

    auto diffuse = Color{0, 0, 0};
    auto specular = Color{0, 0, 0};

    if(lightDotNormal >= 0){
        //Compute the diffuse contribuition
        diffuse = effectiveColor * material.diffuse * lightDotNormal;

        //reflectDotEye represents the cosine of the angle between the reflection vector and the eye vector.
        //A negative number means the light reflects away from the eye.
        auto reflectv = reflect(-lightv, normalv);
        auto reflectDotEye = dot(reflectv, eyev);

        if(reflectDotEye > 0){
            //Compute the specular contribuition
            auto factor = pow(reflectDotEye, material.shininess);
            specular = light.intensity * material.specular * factor;
        }
    }

    return ambient + diffuse + specular;
}

void drawSphereRaycast(unsigned int canvasPixels){

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

void drawSpherePhong(unsigned int canvasPixels){

    auto c = canvas(canvasPixels, canvasPixels);
    auto s = sphere();
    s.material.color = Color{ 1, 0.2, 1 };

    auto light = pointLight(point(-10, 10, -10), Color{1, 1, 1});

    auto rayOrigin = point(0, 0, -5);
    float wallSize = 7;
    float wallZ = 10;
    float pixelSize = wallSize / canvasPixels;
    float half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto ray = Ray{rayOrigin, normalize(point(worldX, worldY, wallZ) - rayOrigin)};
            ray.direction = normalize(ray.direction);

            auto closest = hit(intersect(ray, s));

            if(closest.time >= 0){
                auto point = position(ray, closest.time);
                auto normal = normalAt(s, point);
                auto eye = - ray.direction;

                auto color = lighting(s.material, light, point, eye, normal);
                writePixel(&c, x, y, color);
            }
        }
    }

    canvasToPNG(&c, "phong.png");
}

#endif //RAYTRACERCHALLENGE_RAY_H
