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

struct Bounds {
    Tuple min;
    Tuple max;
};

struct Object {
    unsigned int id;
    Matrix4x4 transform;
    Material material;
    Object* parent = nullptr;
    Bounds bounds;

    Object(){
        static unsigned int cont = 0;
        id = cont++;
        transform = Matrix4x4();
        material = Material();
        bounds = Bounds{point(-INFINITY, -INFINITY, -INFINITY), point(INFINITY, INFINITY, INFINITY)};
    }

    virtual ~Object() = default;
};

Bounds boundsFor(Object* object);

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

struct Cylinder : Object {
    double minimumY;
    double maximumY;
    bool closed;

    Cylinder(): Object() {
        minimumY = 0;
        maximumY = 2;
        closed = false;
    }
};

struct Cone : Cylinder {};

struct Group : Object{
    std::vector<Sphere> spheres = std::vector<Sphere>();
    std::vector<Plane> planes = std::vector<Plane>();
    std::vector<Cube> cubes = std::vector<Cube>();
    std::vector<Cylinder> cylinders = std::vector<Cylinder>();
    std::vector<Cone> cones = std::vector<Cone>();
    std::vector<Group> groups = std::vector<Group>();

    unsigned long long int size(){
        return spheres.size() + planes.size() + cubes.size() + cylinders.size() + cones.size() + groups.size();
    }

    Object* get(unsigned int position){
        if(position < spheres.size()){
            return &spheres[position];
        }
        position -= spheres.size();

        if(position < planes.size()){
            return &planes[position];
        }
        position -= planes.size();

        if(position < cubes.size()){
            return &cubes[position];
        }
        position -= cubes.size();

        if(position < cylinders.size()){
            return &cylinders[position];
        }
        position -= cylinders.size();

        if(position < cones.size()){
            return &cones[position];
        }
        position -= cones.size();

        if(position < groups.size()){
            return &groups[position];
        }

        return nullptr;
    }

    bool insert(Object* object){

        object->parent = this;

        auto pSphere = dynamic_cast<Sphere*>(object);
        if(pSphere != nullptr){
            spheres.push_back(*pSphere);
            return true;
        }

        auto pPlane = dynamic_cast<Plane*>(object);
        if(pPlane != nullptr){
            planes.push_back(*pPlane);
            return true;
        }

        auto pCube = dynamic_cast<Cube*>(object);
        if(pCube != nullptr){
            cubes.push_back(*pCube);
            return true;
        }

        auto pCone = dynamic_cast<Cone*>(object);
        if(pCone != nullptr){
            cones.push_back(*pCone);
            return true;
        }

        auto pCylinder = dynamic_cast<Cylinder*>(object);
        if(pCylinder != nullptr){
            cylinders.push_back(*pCylinder);
            return true;
        }

        auto pGroup = dynamic_cast<Group*>(object);
        if(pGroup != nullptr){
            groups.push_back(*pGroup);
            return true;
        }

        return false;
    }

    void generateBounds(){
        this->bounds = boundsFor(this);
    }
};

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
    std::vector<Object*> objects  = std::vector<Object*>();

    World() = default;
};

struct Tuple2D {
    double x;
    double y;
};

Bounds localBoundsFor(Sphere* sphere){
    return Bounds { point(-1, -1, -1), point(1, 1, 1) };
}

Bounds localBoundsFor(Plane* plane){
    return Bounds { point(-INFINITY, 0, -INFINITY), point(INFINITY, 0, INFINITY) };
}

Bounds localBoundsFor(Cube* cube){
    return Bounds { point(-1, -1, -1), point(1, 1, 1) };
}

Bounds localBoundsFor(Cone* cone){
    return Bounds { point(-1, cone->minimumY, -1), point(1, cone->maximumY, 1) };
}

Bounds localBoundsFor(Cylinder* cylinder){
    return Bounds { point(-1, cylinder->minimumY, -1), point(1, cylinder->maximumY, 1) };
}

Bounds localBoundsFor(Group* group){
    auto bounds = Bounds{ point(INFINITY, INFINITY, INFINITY), point(-INFINITY, -INFINITY, -INFINITY) };

    for(unsigned int i = 0; i < group->size(); i++){
        auto object = group->get(i);
        auto objectBound = boundsFor(object);

        objectBound.min = objectBound.min * object->transform;
        objectBound.max = objectBound.max * object->transform;

        if(lengthSquared(objectBound.min) < lengthSquared(bounds.min)){
            bounds.min = objectBound.min;
        }

        if(lengthSquared(objectBound.max) < lengthSquared(bounds.max)){
            bounds.max = objectBound.max;
        }
    }

    return bounds;
}

Bounds boundsFor(Object* object){
    auto pSphere = dynamic_cast<Sphere*>(object);
    if(pSphere != nullptr){
        return localBoundsFor(pSphere);
    }

    auto pPlane = dynamic_cast<Plane*>(object);
    if(pPlane != nullptr){
        return localBoundsFor(pPlane);
    }

    auto pCube = dynamic_cast<Cube*>(object);
    if(pCube != nullptr){
        return localBoundsFor(pCube);
    }

    auto pCone = dynamic_cast<Cone*>(object);
    if(pCone != nullptr){
        return localBoundsFor(pCone);
    }

    auto pCylinder = dynamic_cast<Cylinder*>(object);
    if(pCylinder != nullptr){
        return localBoundsFor(pCylinder);
    }

    auto pGroup = dynamic_cast<Group*>(object);
    if(pGroup != nullptr){
        return localBoundsFor(pGroup);
    }
}

Tuple2D checkAxis(double origin, double direction, double min, double max){
    auto tMinNumerator = (min - origin);
    auto tMaxNumerator = (max - origin);

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
bool intersectsBounds(Ray ray, Bounds bounds){
    auto xt = checkAxis(ray.origin.x, ray.direction.x, bounds.min.x, bounds.max.x);
    auto yt = checkAxis(ray.origin.y, ray.direction.y, bounds.min.y, bounds.max.y);
    auto zt = checkAxis(ray.origin.z, ray.direction.z, bounds.min.z, bounds.max.z);

    auto tmin = max(xt.x, yt.x, zt.x);
    auto tmax = min(xt.y, yt.y, zt.y);

    return tmin <= tmax;
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

bool checkCap(Ray ray, double t){
    auto x = ray.origin.x + (t * ray.direction.x);
    auto z = ray.origin.z + (t * ray.direction.z);
    return ((x*x) + (z*z)) <= 1;
}

std::vector<Intersection> intersectCaps(Ray ray, Cylinder* cylinder){
    auto lista = std::vector<Intersection>();
    if(!cylinder->closed){
        return lista;
    }

    auto t = (cylinder->minimumY - ray.origin.y) / ray.direction.y;
    if(checkCap(ray, t)){
        lista.push_back(Intersection{cylinder, t});
    }

    t = (cylinder->maximumY - ray.origin.y) / ray.direction.y;
    if(checkCap(ray, t)){
        lista.push_back(Intersection{cylinder, t});
    }

    return lista;
}

std::vector<Intersection> localIntersect(Ray ray, Cylinder* cylinder){
    auto lista = std::vector<Intersection>();

    auto a = (ray.direction.x * ray.direction.x) + (ray.direction.z * ray.direction.z);
    if(!equal(a, 0)){
        auto b = (2 * ray.origin.x * ray.direction.x) + (2 * ray.origin.z * ray.direction.z);
        auto c = (ray.origin.x * ray.origin.x) + (ray.origin.z * ray.origin.z) - 1;

        auto disc = (b*b) - (4 * a * c);

        if(disc < 0){
            return lista;
        }

        auto sqrtDisc = sqrt(disc);
        auto t0 = (-b - sqrtDisc) / (2 * a);
        auto t1 = (-b + sqrtDisc) / (2 * a);

        auto y0 = ray.origin.y + (t0 * ray.direction.y);
        if(cylinder->minimumY < y0 && y0 < cylinder->maximumY){
            lista.push_back(Intersection{ cylinder, t0 });
        }

        auto y1 = ray.origin.y + (t1 * ray.direction.y);
        if(cylinder->minimumY < y1 && y1 < cylinder->maximumY){
            lista.push_back(Intersection{ cylinder, t1 });
        }
    }

    for(auto inter: intersectCaps(ray, cylinder)){
        lista.push_back(inter);
    }

    return lista;
}

bool checkCap(Ray ray, double t, double y){
    auto x = ray.origin.x + (t * ray.direction.x);
    auto z = ray.origin.z + (t * ray.direction.z);
    return ((x*x) + (z*z)) <= y;
}

std::vector<Intersection> intersectCaps(Ray ray, Cone* cylinder){
    auto lista = std::vector<Intersection>();
    if(!cylinder->closed){
        return lista;
    }

    auto t = (cylinder->minimumY - ray.origin.y) / ray.direction.y;
    if(checkCap(ray, t, cylinder->minimumY)){
        lista.push_back(Intersection{cylinder, t});
    }

    t = (cylinder->maximumY - ray.origin.y) / ray.direction.y;
    if(checkCap(ray, t, cylinder->maximumY)){
        lista.push_back(Intersection{cylinder, t});
    }

    return lista;
}

std::vector<Intersection> localIntersect(Ray ray, Cone* cone){
    auto lista = std::vector<Intersection>();

    auto a = (ray.direction.x * ray.direction.x) - (ray.direction.y * ray.direction.y) + (ray.direction.z * ray.direction.z);
    auto b = (2 * ray.origin.x * ray.direction.x) - (2 * ray.origin.y * ray.direction.y) + (2 * ray.origin.z * ray.direction.z);
    auto c = (ray.origin.x * ray.origin.x) - (ray.origin.y * ray.origin.y) + (ray.origin.z * ray.origin.z);

    if(!equal(a, 0)){
        auto disc = (b*b) - (4 * a * c);

        if(disc < 0){
            return lista;
        }

        auto sqrtDisc = sqrt(disc);
        auto t0 = (-b - sqrtDisc) / (2 * a);
        auto t1 = (-b + sqrtDisc) / (2 * a);

        auto y0 = ray.origin.y + (t0 * ray.direction.y);
        if(cone->minimumY < y0 && y0 < cone->maximumY){
            lista.push_back(Intersection{ cone, t0 });
        }

        auto y1 = ray.origin.y + (t1 * ray.direction.y);
        if(cone->minimumY < y1 && y1 < cone->maximumY){
            lista.push_back(Intersection{ cone, t1 });
        }
    }else if(!equal(b, 0)){
        auto t = -c/(2*b);
        lista.push_back(Intersection{ cone, t });
    }

    for(auto inter: intersectCaps(ray, cone)){
        lista.push_back(inter);
    }

    return lista;
}

std::vector<Intersection> intersect(Ray ray, Object* object);

std::vector<Intersection> localIntersect(Ray ray, Group* group){
    auto lista = std::vector<Intersection>();

    if(intersectsBounds(ray, group->bounds)){
        for (unsigned int i = 0; i < group->size(); i++) {
            for(auto intersection: intersect(ray, group->get(i))){
                unsigned int pos = 0;
                for(; pos < lista.size(); pos++){
                    if(lista[pos].time > intersection.time){
                        break;
                    }
                }
                lista.insert(lista.begin() + pos, intersection);
            }
        }
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

    auto pCone = dynamic_cast<Cone*>(object);
    if(pCone != nullptr){
        return localIntersect(localRay, pCone);
    }

    auto pCylinder = dynamic_cast<Cylinder*>(object);
    if(pCylinder != nullptr){
        return localIntersect(localRay, pCylinder);
    }

    auto pGroup = dynamic_cast<Group*>(object);
    if(pGroup != nullptr){
        return localIntersect(localRay, pGroup);
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

Tuple localNormalAt(Cylinder* cylinder, Tuple p){
    auto dist = (p.x*p.x) + (p.z*p.z);

    if(dist < 1){
        if(p.y >= cylinder->maximumY - EPSILON){
            return vector(0, 1, 0);
        }else if(p.y <= cylinder->minimumY + EPSILON){
            return vector(0, -1, 0);
        }
    }

    return vector(p.x, 0, p.z);
}

Tuple localNormalAt(Cone* cone, Tuple p){
    auto dist = (p.x*p.x) + (p.z*p.z);

    if(dist < 1){
        if(p.y >= cone->maximumY - EPSILON){
            return vector(0, 1, 0);
        }else if(p.y <= cone->minimumY + EPSILON){
            return vector(0, -1, 0);
        }
    }

    auto y = sqrt(dist);
    if(p.y > 0){
        y = -y;
    }

    return vector(p.x, y, p.z);
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

    auto pCone = dynamic_cast<Cone*>(object);
    if(pCone != nullptr){
        return localNormalAt(pCone, p);
    }

    auto pCylinder = dynamic_cast<Cylinder*>(object);
    if(pCylinder != nullptr){
        return localNormalAt(pCylinder, p);
    }

    return vector(0, 0, 0);
}

Tuple worldToObject(Object* object, Tuple point){
    if(object->parent != nullptr){
        point = worldToObject(object->parent, point);
    }
    return inverse(object->transform) * point;
}

Tuple normalToWorld(Object* object, Tuple normal){
    normal = transpose(inverse(object->transform)) * normal;
    normal.w = 0;
    normal = normalize(normal);

    if(object->parent != nullptr){
        normal = normalToWorld(object->parent, normal);
    }

    return normal;
}

Tuple normalAt(Object* object, Tuple p){
    auto localPoint = worldToObject(object, p);
    auto localNormal = localNormalAt(object, localPoint);
    return normalToWorld(object, localNormal);
}

#endif //RAYTRACERCHALLENGE_STRUCTS_H
