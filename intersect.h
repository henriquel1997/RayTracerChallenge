//
// Created by Henrique on 22/06/2019.
//

#ifndef RAYTRACERCHALLENGE_INTERSECT_H
#define RAYTRACERCHALLENGE_INTERSECT_H

#include <vector>
#include "structs.h"
#include "ray.h"
#include "bounding_box.h"

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

std::vector<Intersection> localIntersect(Ray ray, Triangle* triangle){
    auto lista = std::vector<Intersection>();

    auto dirCrossE2 = cross(ray.direction, triangle->e2);
    auto det = dot(triangle->e1, dirCrossE2);

    if(absolute(det) >= EPSILON){
        auto f = 1.0 / det;
        auto p1ToOrigin = ray.origin - triangle->p1;
        auto u = f * dot(p1ToOrigin, dirCrossE2);

        if(u >= 0 && u <= 1){
            auto originCrossE1 = cross(p1ToOrigin, triangle->e1);
            auto v = f * dot(ray.direction, originCrossE1);

            if(v >= 0 && (u + v) <= 1){
                auto t = f * dot(triangle->e2, originCrossE1);
                lista.push_back(Intersection{triangle, t, u, v});
            }
        }
    }

    return lista;
}

bool intersectionAllowed(OperationCSG operation, bool leftHit, bool insideLeft, bool insideRight){
    if(operation == UNION){
        return (leftHit && !insideRight) || (!leftHit && !insideLeft);
    }else if(operation == INTERSECTION){
        return (leftHit && insideRight) || (!leftHit && insideLeft);
    }else if(operation == DIFFERENCE){
        return (leftHit && !insideRight) || (!leftHit && insideLeft);
    }
    return false;
}

std::vector<Intersection> filterIntersections(CSG* csg, const std::vector<Intersection> &intersections){
    bool insideLeft = false;
    bool insideRight = false;

    auto result = std::vector<Intersection>();

    for(auto i : intersections){
        auto leftHit = csg->left->includes(i.object);

        if(intersectionAllowed(csg->operation, leftHit, insideLeft, insideRight)){
            result.push_back(i);
        }

        if(leftHit){
            insideLeft = !insideLeft;
        }else{
            insideRight = !insideRight;
        }
    }

    return result;
}

std::vector<Intersection> intersect(Ray ray, Object* object);

std::vector<Intersection> localIntersect(Ray ray, CSG* csg){
    auto result = std::vector<Intersection>();

    //Chca se o raio intercede o CSG
    if(intersects(ray, boundsOf(csg))){
        //Insere as interseções do lado esquerdo em ordem
        for(auto intersection : intersect(ray, csg->left)){
            unsigned int pos = 0;
            for(; pos < result.size(); pos++){
                if(result[pos].time > intersection.time){
                    break;
                }
            }
            result.insert(result.begin() + pos, intersection);
        }

        //Insere as interseções do lado direito em ordem
        for(auto intersection : intersect(ray, csg->right)){
            unsigned int pos = 0;
            for(; pos < result.size(); pos++){
                if(result[pos].time > intersection.time){
                    break;
                }
            }
            result.insert(result.begin() + pos, intersection);
        }

        //Retorna o resultado filtrado
        return filterIntersections(csg, result);
    }

    return result;
}

std::vector<Intersection> localIntersect(Ray ray, Group* group){
    auto lista = std::vector<Intersection>();
    if(intersects(ray, boundsOf(group))){
        for (unsigned long long i = 0; i < group->size(); i++) {
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

    auto pTriangle = dynamic_cast<Triangle*>(object);
    if(pTriangle != nullptr){
        return localIntersect(localRay, pTriangle);
    }

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

    auto pCSG = dynamic_cast<CSG*>(object);
    if(pCSG != nullptr){
        return localIntersect(localRay, pCSG);
    }

    return std::vector<Intersection>();
}

#endif //RAYTRACERCHALLENGE_INTERSECT_H
