//
// Created by Henrique on 07/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RENDERING_H
#define RAYTRACERCHALLENGE_RENDERING_H

#include "matrix.h"
#include "color.h"
#include "canvas.h"
#include <vector>
#include "ray.h"
#include "transformations.h"
#include "pattern.h"
#include "structs.h"

struct Sphere : Object {
    Tuple origin{};
    double radius;

    Sphere(): Object(){
        origin = point(0, 0, 0);
        radius = 1;
    }
};

struct Plane : Object {};

struct Intersection{
    Object* object;
    double time;
};

struct Computations{
    Intersection intersection;
    Tuple point;
    Tuple eyev;
    Tuple normalv;
    bool inside;
    Tuple overPoint;
};

struct Camera{
    unsigned int hSize;
    unsigned int vSize;
    double fieldOfView;
    double pixelSize;
    double halfWidth;
    double halfHeight;
    Matrix4x4 transform;

    Camera(unsigned int hSize, unsigned int vSize, double fieldOfView){
        this->hSize = hSize;
        this->vSize = vSize;
        this->fieldOfView = fieldOfView;
        this->transform = Matrix4x4();

        auto halfView = tan(fieldOfView/2);
        auto aspectRatio = hSize / vSize;

        if(aspectRatio >= 1){
            this->halfWidth = halfView;
            this->halfHeight = halfView/aspectRatio;
        }else{
            this->halfWidth = halfView/aspectRatio;
            this->halfHeight = halfView;
        }

        this->pixelSize = (this->halfWidth * 2) / hSize;
    }
};

struct World{
    std::vector<Light> lightSources = std::vector<Light>();
    std::vector<Object*> objects  = std::vector<Object*>();;

    World() = default;
};

void addPattern(Material* material, Pattern* pattern){
    material->hasPattern = true;
    material->pattern = pattern;
}

//Expects a local ray (needs to go through transform(ray, inverse(sphere.transform)))
std::vector<Intersection> localIntersect(Ray ray, Sphere* sphere){
    auto sphere_to_ray = ray.origin - sphere->origin;

    auto a = dot(ray.direction, ray.direction);
    auto b = 2 * dot(ray.direction, sphere_to_ray);
    auto c = dot(sphere_to_ray, sphere_to_ray) - 1;

    auto discriminant = (b * b) - 4 * a * c;

    if(discriminant < 0){
        return std::vector<Intersection>();
    }

    auto raiz = sqrt(discriminant);
    auto divisor = (2 * a);

    auto t1 = (-b - raiz) / divisor;
    auto t2 = (-b + raiz) / divisor;

    auto lista = std::vector<Intersection>();

    lista.push_back(Intersection{sphere, t1});
    lista.push_back(Intersection{sphere, t2});

    return lista;
}

std::vector<Intersection> localIntersect(Ray ray, Plane* plane){
    auto lista = std::vector<Intersection>();
    if(absolute(ray.direction.y) >= EPSILON){
        auto time = (-ray.origin.y) / ray.direction.y;
        lista.push_back(Intersection{plane, time});
    }
    return lista;
}

std::vector<Intersection> intersect(Ray ray, Object* object){
    auto localRay = transform(ray, inverse(object->transform));

    auto pSphere = dynamic_cast<Sphere*>(object);
    if(pSphere != nullptr){
        return localIntersect(localRay, pSphere);
    }

    auto pPlane = dynamic_cast<Plane*>(object);
    if(pPlane != nullptr){
        return localIntersect(localRay, pPlane);
    }

    return std::vector<Intersection>();
}

Intersection hit(const std::vector<Intersection>& intersections){
    auto hit = Intersection{nullptr, -1};
    for(const auto &intersection: intersections){
        if(intersection.time >= 0 && (intersection.time < hit.time || hit.time == -1)){
            hit = intersection;
        }
    }
    return hit;
}

//Expects a local point (needs to go through inverse(sphere->transform) * point)
Tuple localNormalAt(Sphere* sphere, Tuple p){
    return p - sphere->origin;
}

Tuple localNormalAt(Plane* plane){
    return vector(0, 1, 0);
}

Tuple normalAt(Object* object, Tuple p){
    auto localPoint = inverse(object->transform) * p;
    auto localNormal = vector(0, 0, 0);

    auto pSphere = dynamic_cast<Sphere*>(object);
    if(pSphere != nullptr){
        localNormal = localNormalAt(pSphere, localPoint);
    }else{
        auto pPlane = dynamic_cast<Plane*>(object);
        if(pPlane != nullptr){
            localNormal = localNormalAt(pPlane);
        }
    }

    auto worldNormal = transpose(inverse(object->transform)) * localNormal;
    worldNormal.w = 0;

    return normalize(worldNormal);
}

Color lighting(Material material, Object* object, Light light, Tuple point, Tuple eyev, Tuple normalv, bool inShadow){

    Color color{};
    if(material.hasPattern){
        color = patternAtObject(material.pattern, object, point);
    }else{
        color = material.color;
    }

    //Combine the surface color with the light's color/intensity
    auto effectiveColor = color * light.intensity;

    //Compute the ambient contribution
    auto ambient = effectiveColor * material.ambient;

    if(inShadow){
        return ambient;
    }

    //Find the direction to the light source
    auto lightv = normalize(light.position - point);

    //lightDotNormal represents the cosine of the angle between the light vector and the normal vector.
    //A negative number means the light is on the other side of the surface.
    auto lightDotNormal = dot(lightv, normalv);

    auto diffuse = BLACK;
    auto specular = BLACK;

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

World defaultWorld(){
    auto world = World();

    auto color = WHITE;
    auto defaulLight = pointLight(Tuple{-10, 10, -10}, color);
    world.lightSources.push_back(defaulLight);

    auto s1 = Sphere();
    s1.material.color = Color{0.8, 1.0, 0.6};
    s1.material.diffuse = 0.7;
    s1.material.specular = 0.2;
    world.objects.push_back(&s1);

    auto s2 = Sphere();
    s2.transform = scaling(0.5, 0.5, 0.5);
    world.objects.push_back(&s2);

    return world;
}

std::vector<Intersection> intersectWorld(Ray ray, World world){
    auto lista = std::vector<Intersection>();
    for(Object* object: world.objects){
        auto intersecoes = intersect(ray, object);
        for(const auto &intersecao: intersecoes){
            unsigned int pos = 0;
            for(; pos < lista.size(); pos++){
                if(lista[pos].time > intersecao.time){
                    break;
                }
            }
            lista.insert(lista.begin() + pos, intersecao);
        }
    }
    return lista;
}

Computations prepareComputations(Ray ray, Intersection intersection){
    auto comps = Computations{};
    comps.intersection = intersection;
    comps.point = position(ray, intersection.time);
    comps.eyev = - ray.direction;
    comps.normalv = normalAt(intersection.object, comps.point);
    comps.inside = dot(comps.normalv, comps.eyev) < 0;

    if(comps.inside){
        comps.normalv = - comps.normalv;
    }

    comps.overPoint = comps.point + (comps.normalv * EPSILON);

    return comps;
}

bool isShadowed(World world, Tuple point){
    for(auto light: world.lightSources){
        auto v = light.position - point;
        auto distance = length(point);
        auto direction = normalize(v);

        auto r = Ray{ point, direction };
        auto intersections = intersectWorld(r, world);
        auto h = hit(intersections);
        if(h.time >= 0 && h.time < distance){
            return true;
        }
    }
    return false;
}

Color shadeHit(World world, Computations comps, bool shadows){
    auto material = comps.intersection.object->material;

    auto color = BLACK;
    for(auto light: world.lightSources){
        auto shadowed = shadows && isShadowed(world, comps.overPoint);
        color = color + lighting(material, comps.intersection.object, light, comps.overPoint, comps.eyev, comps.normalv, shadowed);
    }
    return color;
}

Color colorAt(Ray ray, const World& world, bool shadows){
    auto intersections = intersectWorld(ray, world);
    auto h = hit(intersections);
    auto color = BLACK;

    if(h.time >= 0){
        auto comps = prepareComputations(ray, h);
        color = shadeHit(world, comps, shadows);
    }

    return color;
}

Matrix4x4 viewTransform(Tuple from, Tuple to, Tuple up){
    auto forward = normalize(to - from);
    auto left = cross(forward, normalize(up));
    auto trueUp = cross(left, forward);
    auto orientation = Matrix4x4(left.x,     left.y,     left.z,     0,
                                 trueUp.x,   trueUp.y,   trueUp.z,   0,
                                 -forward.x, -forward.y, -forward.z, 0,
                                 0,          0,          0,          1);
    return orientation * translation(-from.x, -from.y, -from.z);
}

Ray rayForPixel(Camera camera, unsigned int x, unsigned int y){
    auto xOffset = (x + 0.5f) * camera.pixelSize;
    auto yOffset = (y + 0.5f) * camera.pixelSize;

    auto worldX = camera.halfWidth - xOffset;
    auto worldY = camera.halfHeight - yOffset;

    auto pixel = inverse(camera.transform) * point(worldX, worldY, -1);
    auto origin = inverse(camera.transform) * point(0, 0, 0);
    auto direction = normalize(pixel - origin);

    return Ray{origin, direction};
}

Canvas render(Camera camera, const World& world, bool shadows){
    auto image = canvas(camera.hSize, camera.vSize);

    for(unsigned int y = 0; y < camera.vSize; y++){
        for(unsigned int x = 0; x < camera.hSize; x++){
            auto ray = rayForPixel(camera, x, y);
            auto color = colorAt(ray, world, shadows);
            writePixel(&image, x, y, color);
        }
    }

    return image;
}

#endif //RAYTRACERCHALLENGE_RENDERING_H
