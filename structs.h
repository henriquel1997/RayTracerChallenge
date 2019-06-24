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
    Object* parent = nullptr;

    Object(){
        static unsigned int cont = 0;
        id = cont++;
        transform = Matrix4x4();
        material = Material();
    }

    virtual ~Object() = default;

    bool includes(Object* o);
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

struct Triangle : Object{
    Tuple p1{};
    Tuple p2{};
    Tuple p3{};

    Tuple e1{};
    Tuple e2{};

    bool smooth;

    Tuple n1{};
    Tuple n2{};
    Tuple n3{};

    Triangle(Tuple p1, Tuple p2, Tuple p3): Object() {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        this->e1 = p2 - p1;
        this->e2 = p3 - p1;
        this->n1 = normalize(cross(this->e2, this->e1));
        this->n2 = this->n1;
        this->n3 = this->n1;
        this->smooth = false;
    }

    Triangle(Tuple p1, Tuple p2, Tuple p3, Tuple n1, Tuple n2, Tuple n3): Object() {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        this->e1 = p2 - p1;
        this->e2 = p3 - p1;
        this->n1 = n1;
        this->n2 = n2;
        this->n3 = n3;
        this->smooth = true;
    }
};

enum OperationCSG {
    UNION, INTERSECTION, DIFFERENCE
};

struct CSG : Object {
    OperationCSG operation;
    Object* left;
    Object* right;

    CSG(OperationCSG operation, Object* left, Object* right) : Object(){
        this->operation = operation;
        left->parent = this;
        right->parent = this;
        this->left = left;
        this->right = right;
    }
};

struct Group : Object{
    std::vector<Triangle> triangles = std::vector<Triangle>();
    std::vector<Sphere> spheres = std::vector<Sphere>();
    std::vector<Plane> planes = std::vector<Plane>();
    std::vector<Cube> cubes = std::vector<Cube>();
    std::vector<Cylinder> cylinders = std::vector<Cylinder>();
    std::vector<Cone> cones = std::vector<Cone>();
    std::vector<Group> groups = std::vector<Group>();
    std::vector<CSG> csgs = std::vector<CSG>();

    unsigned long long int size(){
        return triangles.size() + spheres.size() + planes.size() + cubes.size() + cylinders.size() + cones.size() + groups.size();
    }

    Object* get(unsigned long long position){

        if(position < triangles.size()){
            return &triangles[position];
        }
        position -= triangles.size();

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
        position -= groups.size();

        if(position < csgs.size()){
            return &csgs[position];
        }

        return nullptr;
    }

    bool insert(Object* object){

        object->parent = this;

        auto pTriangle = dynamic_cast<Triangle*>(object);
        if(pTriangle != nullptr){
            triangles.push_back(*pTriangle);
            return true;
        }

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

        auto pCSG = dynamic_cast<CSG*>(object);
        if(pCSG != nullptr){
            csgs.push_back(*pCSG);
            return true;
        }

        return false;
    }

    void remove(unsigned long long position){

        if(position < triangles.size()){
            triangles.erase(triangles.begin() + position);
            return;
        }
        position -= triangles.size();

        if(position < spheres.size()){
            spheres.erase(spheres.begin() + position);
            return;
        }
        position -= spheres.size();

        if(position < planes.size()){
            planes.erase(planes.begin() + position);
            return;
        }
        position -= planes.size();

        if(position < cubes.size()){
            cubes.erase(cubes.begin() + position);
            return;
        }
        position -= cubes.size();

        if(position < cylinders.size()){
            cylinders.erase(cylinders.begin() + position);
            return;
        }
        position -= cylinders.size();

        if(position < cones.size()){
            cones.erase(cones.begin() + position);
            return;
        }
        position -= cones.size();

        if(position < groups.size()){
            groups.erase(groups.begin() + position);
            return;
        }
        position -= groups.size();

        if(position < csgs.size()){
            csgs.erase(csgs.begin() + position);
            return;
        }
    }
};

bool Object::includes(Object* o){
    auto pGroup = dynamic_cast<Group*>(o);
    if(pGroup != nullptr){
        for(unsigned int i = 0; i < pGroup->size(); i++){
            if(pGroup->get(i)->includes(o)){
                return true;
            }
        }
        return false;
    }

    auto pCSG = dynamic_cast<CSG*>(o);
    if(pCSG != nullptr){
        return pCSG->left->includes(o) || pCSG->right->includes(o);
    }

    //Caso não for Grupo ou CSG
    return this->id == o->id;
}

struct Intersection{
    Object* object;
    double time;
    double u;
    double v;
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



#endif //RAYTRACERCHALLENGE_STRUCTS_H
