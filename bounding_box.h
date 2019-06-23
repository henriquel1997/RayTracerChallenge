//
// Created by Henrique on 22/06/2019.
//

#ifndef RAYTRACERCHALLENGE_BOUNDING_BOX_H
#define RAYTRACERCHALLENGE_BOUNDING_BOX_H

#include "tuples.h"
#include "structs.h"
#include "ray.h"

struct BoundingBox{
    Tuple min{};
    Tuple max{};

    BoundingBox(){
        this->min = point(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
        this->max = point(MINDOUBLE, MINDOUBLE, MINDOUBLE);
    }

    BoundingBox(Tuple min, Tuple max){
        this->min = min;
        this->max = max;
    }

    void add(Tuple point){
        if(point.x < min.x){
            min.x = point.x;
        }

        if(point.x > max.x){
            max.x = point.x;
        }

        if(point.y < min.y){
            min.y = point.y;
        }

        if(point.y > max.y){
            max.y = point.y;
        }

        if(point.z < min.z){
            min.z = point.z;
        }

        if(point.z > max.z){
            max.z = point.z;
        }
    }

    void add(BoundingBox boundingBox){
        this->add(boundingBox.min);
        this->add(boundingBox.max);
    }

    bool containsPoint(Tuple p){
        return p.x >= min.x && p.y >= min.y && p.z >= min.z &&
               p.x <= max.x && p.y <= max.y && p.z <= max.z;
    }

    bool containsBox(BoundingBox* boundingBox){
        bool containsMin = containsPoint(boundingBox->min);
        bool containsMax = containsPoint(boundingBox->max);
        return containsMin && containsMax;
    }

    void transform(Matrix4x4 matrix){
        Tuple points[8];
        points[0] = min;
        points[1] = point(min.x, min.y, max.z);
        points[2] = point(min.x, max.y, min.z);
        points[3] = point(min.x, max.y, max.z);
        points[4] = point(max.x, min.y, min.z);
        points[5] = point(max.x, min.y, max.z);
        points[6] = point(max.x, max.y, min.z);
        points[7] = max;

        this->min = point(MAXDOUBLE, MAXDOUBLE, MAXDOUBLE);
        this->max = point(MINDOUBLE, MINDOUBLE, MINDOUBLE);

        for(auto point: points){
            add(matrix * point);
        }
    }
};

BoundingBox boundsOf(Sphere* sphere){
    auto min = point(-sphere->radius, -sphere->radius, -sphere->radius);
    auto max = point(sphere->radius, sphere->radius, sphere->radius);
    return {min, max};
}

BoundingBox boundsOf(Plane* plane){
    return {point(MINDOUBLE, 0, MINDOUBLE), point(MAXDOUBLE, 0, MAXDOUBLE)};
}

BoundingBox boundsOf(Cube* cube){
    return {point(-1, -1, -1), point(1, 1, 1)};
}

BoundingBox boundsOf(Cylinder* cylinder){
    return {point(-1, cylinder->minimumY, -1), point(1, cylinder->maximumY, 1)};
}

BoundingBox boundsOf(Cone* cone){
    return {point(-1, cone->minimumY, -1), point(1, cone->maximumY, 1)};
}

BoundingBox boundsOf(Triangle* triangle){
    auto bbox = BoundingBox();
    bbox.add(triangle->p1);
    bbox.add(triangle->p2);
    bbox.add(triangle->p3);
    return bbox;
}

BoundingBox boundsOf(Group* group);

BoundingBox boundsOf(Object* object){

    auto pTriangle = dynamic_cast<Triangle*>(object);
    if(pTriangle != nullptr){
        return boundsOf(pTriangle);
    }

    auto pSphere = dynamic_cast<Sphere*>(object);
    if(pSphere != nullptr){
        return boundsOf(pSphere);
    }

    auto pPlane = dynamic_cast<Plane*>(object);
    if(pPlane != nullptr){
        return boundsOf(pPlane);
    }

    auto pCube = dynamic_cast<Cube*>(object);
    if(pCube != nullptr){
        return boundsOf(pCube);
    }

    auto pCone = dynamic_cast<Cone*>(object);
    if(pCone != nullptr){
        return boundsOf(pCone);
    }

    auto pCylinder = dynamic_cast<Cylinder*>(object);
    if(pCylinder != nullptr){
        return boundsOf(pCylinder);
    }

    auto pGroup = dynamic_cast<Group*>(object);
    if(pGroup != nullptr){
        return boundsOf(pGroup);
    }

    return BoundingBox();
}

BoundingBox parentSpaceBoundsOf(Object* object){
    auto bbox = boundsOf(object);
    bbox.transform(object->transform);
    return bbox;
}

BoundingBox boundsOf(Group* group){
    auto bbox = BoundingBox();
    for(unsigned int i = 0; i < group->size(); i++){
        bbox.add(parentSpaceBoundsOf(group->get(i)));
    }
    return bbox;
}

struct Point2D {
    double x;
    double y;
};

Point2D checkAxis(double origin, double direction, double min, double max){
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

    return Point2D{ tmin, tmax };
}

bool intersects(Ray ray, BoundingBox boundingBox){
    auto xt = checkAxis(ray.origin.x, ray.direction.x, boundingBox.min.x, boundingBox.max.x);
    auto yt = checkAxis(ray.origin.y, ray.direction.y, boundingBox.min.y, boundingBox.max.y);
    auto zt = checkAxis(ray.origin.z, ray.direction.z, boundingBox.min.z, boundingBox.max.z);

    auto tmin = max(xt.x, yt.x, zt.x);
    auto tmax = min(xt.y, yt.y, zt.y);

    return tmin <= tmax;
}

struct Halfs{
    BoundingBox half1;
    BoundingBox half2;

};

Halfs splitBounds(BoundingBox boundingBox){
    auto dx = boundingBox.max.x - boundingBox.min.x;
    auto dy = boundingBox.max.y - boundingBox.min.y;
    auto dz = boundingBox.max.z - boundingBox.min.z;

    auto greatest = max(dx, dy, dz);

    auto x0 = boundingBox.min.x;
    auto y0 = boundingBox.min.y;
    auto z0 = boundingBox.min.z;
    auto x1 = boundingBox.max.x;
    auto y1 = boundingBox.max.y;
    auto z1 = boundingBox.max.z;

    if(greatest == dx){
        x0 = x0 + (dx/2);
        x1 = x0;
    }else if(greatest == dy){
        y0 = y0 + (dx/2);
        y1 = y0;
    }else{
        z0 = z0 + (dx/2);
        z1 = z0;
    }

    auto midMin = point(x0, y0, z0);
    auto midMax = point(x1, y1, z1);

    auto left = BoundingBox(boundingBox.min, midMax);
    auto right = BoundingBox(midMin, boundingBox.max);

    return Halfs{left, right};
}

void makeSubGroups(Group* group){
    auto halfs = splitBounds(parentSpaceBoundsOf(group));
    auto group1 = new Group();
    auto group2 = new Group();

    for(unsigned int i = 0; i < group->size(); i++){
        auto object = group->get(i);
        auto bound = parentSpaceBoundsOf(object);

        if(halfs.half1.containsBox(&bound)){
            group1->insert(object);
            group->remove(i);
            i--;
        }else if(halfs.half2.containsBox(&bound)){
            group2->insert(object);
            group->remove(i);
            i--;
        }
    }

    if(group1->size() > 0){
        group->insert(group1);
    }

    if(group2->size() > 0){
        group->insert(group2);
    }
}

void divide(Group* group, unsigned int threshold){
    if(group->size() >= threshold){
        makeSubGroups(group);
    }

    for (auto &subGroup : group->groups) {
        divide(&subGroup, threshold);
    }
}

#endif //RAYTRACERCHALLENGE_BOUNDING_BOX_H
