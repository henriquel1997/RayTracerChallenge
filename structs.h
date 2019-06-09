//
// Created by Henrique on 08/06/2019.
//

#ifndef RAYTRACERCHALLENGE_STRUCTS_H
#define RAYTRACERCHALLENGE_STRUCTS_H


#include "matrix.h"
#include "color.h"

struct Pattern{
    Color a{};
    Color b{};
    Matrix4x4 transform;

    Pattern() = default;

    Pattern(Color a, Color b){
        this->a = a;
        this->b = b;
        transform = identity();
    }
    virtual ~Pattern() = default;
};

struct Material{
    Color color{};
    double ambient;
    double diffuse;
    double specular;
    double shininess;
    bool hasPattern;
    Pattern* pattern;
    double reflective;
    double transparency;
    double refractiveIndex;
    bool castShadows;

    Material(){
        color = WHITE;
        ambient = 0.1f;
        diffuse = 0.9f;
        specular = 0.9f;
        shininess = 200;
        hasPattern = false;
        pattern = nullptr;
        reflective = 0;
        transparency = 0;
        refractiveIndex = 1;
        castShadows = true;
    }
};

struct Object {
    unsigned int id;
    Matrix4x4 transform;
    Material material;

    Object(){
        static unsigned int cont = 0;
        id = cont++;
        transform = Matrix4x4();
        material = Material();
    }

    virtual ~Object() = default;
};

struct Sphere : Object {
    Tuple origin{};
    double radius;

    Sphere(): Object(){
        origin = point(0, 0, 0);
        radius = 1;
    }
};

struct Plane : Object {};

struct Cube : Object {};

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
    Tuple underPoint;
    Tuple reflectv;
    double n1; //For Refraction, refractive index of the material the ray is passing from
    double n2; //For Refraction, refractive index of the material the ray is passing to
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

struct Tuple2D {
    double x;
    double y;
};

Tuple2D checkAxis(double origin, double direction){
    auto tMinNumerator = (-1 - origin);
    auto tMaxNumerator = (1 - origin);

    double tmin, tmax;
    if(absolute(direction) >= EPSILON){
        tmin = tMinNumerator / direction;
        tmax = tMaxNumerator / direction;
    }else{
        tmin = tMinNumerator * INFINITY;
        tmax = tMaxNumerator * INFINITY;
    }

    if(tmin > tmax){
        double aux = tmin;
        tmin = tmax;
        tmax = aux;
    }

    return Tuple2D{ tmin, tmax };
}

//TODO: Pode ser otimizado
std::vector<Intersection> localIntersect(Ray ray, Cube* cube){
    auto xt = checkAxis(ray.origin.x, ray.direction.x);
    auto yt = checkAxis(ray.origin.y, ray.direction.y);
    auto zt = checkAxis(ray.origin.z, ray.direction.z);

    auto tmin = max(xt.x, yt.x, zt.x);
    auto tmax = min(xt.y, yt.y, zt.y);

    auto lista = std::vector<Intersection>();

    if(tmin <= tmax){
        lista.push_back(Intersection{ cube, tmin });
        lista.push_back(Intersection{ cube, tmax });
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

    auto pCube = dynamic_cast<Cube*>(object);
    if(pCube != nullptr){
        return localIntersect(localRay, pCube);
    }

    return std::vector<Intersection>();
}

//Expects a local point (needs to go through inverse(sphere->transform) * point)
Tuple localNormalAt(Sphere* sphere, Tuple p){
    return p - sphere->origin;
}

Tuple localNormalAt(Plane* plane){
    return vector(0, 1, 0);
}

Tuple localNormalAt(Cube* cube, Tuple p){
    auto absX = absolute(p.x);
    auto absY = absolute(p.y);
    auto absZ = absolute(p.z);

    auto maxc = max(absX, absY, absZ);

    if(maxc == absX){
        return vector(p.x, 0, 0);
    }else if(maxc == absY){
        return vector(0, p.y, 0);
    }
    return vector(0, 0, p.z);
}

Tuple localNormalAt(Object* object, Tuple p){
    auto pSphere = dynamic_cast<Sphere*>(object);
    if(pSphere != nullptr){
        return localNormalAt(pSphere, p);
    }

    auto pPlane = dynamic_cast<Plane*>(object);
    if(pPlane != nullptr){
        return localNormalAt(pPlane);
    }

    auto pCube = dynamic_cast<Cube*>(object);
    if(pCube != nullptr){
        return localNormalAt(pCube, p);
    }

    return vector(0, 0, 0);
}

Tuple normalAt(Object* object, Tuple p){
    auto localPoint = inverse(object->transform) * p;
    auto worldNormal = transpose(inverse(object->transform)) * localNormalAt(object, localPoint);
    worldNormal.w = 0;
    return normalize(worldNormal);
}


#endif //RAYTRACERCHALLENGE_STRUCTS_H
