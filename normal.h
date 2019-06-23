//
// Created by Henrique on 22/06/2019.
//

#ifndef RAYTRACERCHALLENGE_NORMAL_H
#define RAYTRACERCHALLENGE_NORMAL_H

#include "structs.h"

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

Tuple localNormalAt(Triangle* triangle, Intersection hit){
    if(triangle->smooth){
        return (triangle->n2 * hit.u) + (triangle->n3 * hit.v) + (triangle->n1 * (1 - hit.u - hit.v));
    }
    return triangle->n1;
}

Tuple localNormalAt(Object* object, Tuple p, Intersection hit){

    auto pTriangle = dynamic_cast<Triangle*>(object);
    if(pTriangle != nullptr){
        return localNormalAt(pTriangle, hit);
    }

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

Tuple normalAt(Object* object, Tuple p, Intersection hit){
    auto localPoint = worldToObject(object, p);
    auto localNormal = localNormalAt(object, localPoint, hit);
    return normalToWorld(object, localNormal);
}

#endif //RAYTRACERCHALLENGE_NORMAL_H
