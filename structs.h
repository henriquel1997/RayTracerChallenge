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

    Material(){
        color = WHITE;
        ambient = 0.1f;
        diffuse = 0.9f;
        specular = 0.9f;
        shininess = 200;
        hasPattern = false;
        pattern = nullptr;
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
};;

#endif //RAYTRACERCHALLENGE_STRUCTS_H
