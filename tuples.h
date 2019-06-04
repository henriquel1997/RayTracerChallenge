//
// Created by Henrique on 04/06/2019.
//

#ifndef RAYTRACERCHALLENGE_TUPLES_H
#define RAYTRACERCHALLENGE_TUPLES_H

#include <cstdio>
#include <cmath>
#include "util.h"

struct Tuple{
    float x;
    float y;
    float z;
    float w;
};

Tuple point(float x, float y, float z){
    return Tuple{ x, y, z, 1.f};
}

Tuple vector(float x, float y, float z){
    return Tuple{ x, y, z, 0.f};
}

bool isPoint(Tuple* t){
    return t->w == 1.f;
}

bool isVector(Tuple* t){
    return t->w == 0.f;
}

bool equal(Tuple* t1, Tuple* t2){
    return equal(t1->x, t2->x) &&
           equal(t1->y, t2->y) &&
           equal(t1->z, t2->z) &&
           equal(t1->w, t2->w);
}

bool operator == (Tuple t1, Tuple t2) {
    return equal(&t1, &t2);
}

Tuple add(Tuple* t1, Tuple* t2){
    return Tuple{ t1->x + t2->x,
                  t1->y + t2->y,
                  t1->z + t2->z,
                  t1->w + t2->w };
}

Tuple operator + (Tuple t1, Tuple t2){
    return add(&t1, &t2);
}

Tuple subtract(Tuple* t1, Tuple* t2){
    return Tuple{ t1->x - t2->x,
                  t1->y - t2->y,
                  t1->z - t2->z,
                  t1->w - t2->w };
}

Tuple operator - (Tuple t1, Tuple t2){
    return subtract(&t1, &t2);
}

Tuple negate(Tuple* t1){
    return Tuple{ -t1->x,
                  -t1->y,
                  -t1->z,
                  -t1->w };
}

Tuple operator - (Tuple t){
    return negate(&t);
}

Tuple multiplication(Tuple* t, float f){
    return Tuple{ t->x * f,
                  t->y * f,
                  t->z * f,
                  t->w * f };
}

Tuple operator * (Tuple t, float f){
    return multiplication(&t, f);
}

Tuple operator * (float f, Tuple t){
    return multiplication(&t, f);
}

Tuple division(Tuple* t, float f){
    return Tuple{ t->x / f,
                  t->y / f,
                  t->z / f,
                  t->w / f };
}

Tuple operator / (Tuple t, float f){
    return division(&t, f);
}

float length(Tuple* t){
    return sqrtf((t->x * t->x) + (t->y * t->y) + (t->z * t->z) + (t->w * t->w));
}

float lengthSquared(Tuple* t){
    return (t->x * t->x) + (t->y * t->y) + (t->z * t->z) + (t->w * t->w);
}

Tuple normalize(Tuple* t){
    float magnitude = length(t);
    return Tuple{ t->x / magnitude,
                  t->y / magnitude,
                  t->z / magnitude,
                  t->w / magnitude};
}

float dot(Tuple* t1, Tuple* t2){
    return (t1->x * t2->x) + (t1->y * t2->y) + (t1->z * t2->z) + (t1->w * t2->w);
}

Tuple cross(Tuple* t1, Tuple* t2){
    return vector((t1->y * t2->z) - (t1->z * t2->y),
                  (t1->z * t2->x) - (t1->x * t2->z),
                  (t1->x * t2->y) - (t1->y * t2->x));
}

void printTuple(Tuple* tuple){
    printf("x: %f, y: %f, z: %f, w: %f", tuple->x, tuple->y, tuple->z, tuple->w);
}

#endif //RAYTRACERCHALLENGE_TUPLES_H
