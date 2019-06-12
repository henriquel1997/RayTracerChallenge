//
// Created by Henrique on 08/06/2019.
//

#ifndef RAYTRACERCHALLENGE_PATTERN_H
#define RAYTRACERCHALLENGE_PATTERN_H

#include "color.h"
#include "structs.h"

//TODO: Radial Gradient, Perlin Noise
struct Stripes: Pattern{ Stripes(Color a, Color b): Pattern(a, b){} };
struct Gradient: Pattern{ Gradient(Color a, Color b): Pattern(a, b){} };
struct Ring: Pattern{ Ring(Color a, Color b): Pattern(a, b){} };
struct Checker3D: Pattern{ Checker3D(Color a, Color b): Pattern(a, b){} };
struct NestedChecker3D: Pattern{
    Pattern a;
    Pattern b;
    bool useColors;

    NestedChecker3D(Color a, Color b, const Pattern &p1, const Pattern &p2): Pattern(a, b){
        this->a = p1;
        this->b = p2;
        useColors = true;
    }

    NestedChecker3D(const Pattern &p1, const Pattern &p2): Pattern(){
        this->a = p1;
        this->b = p2;
        useColors = false;
    }
};

struct Blend: Pattern{
    Pattern a;
    Pattern b;
    double amount;
    bool useColors;

    Blend(Color a, Color b, const Pattern &p1, const Pattern &p2): Pattern(a, b){
        this->a = p1;
        this->b = p2;
        useColors = true;
        amount = 0.5;
    }

    Blend(Color a, Color b, const Pattern &p1, const Pattern &p2, double amount): Pattern(a, b){
        this->a = p1;
        this->b = p2;
        useColors = true;
        this->amount = amount;
    }

    Blend(const Pattern &p1, const Pattern &p2): Pattern(){
        this->a = p1;
        this->b = p2;
        useColors = false;
        amount = 0.5;
    }

    Blend(const Pattern &p1, const Pattern &p2, double amount): Pattern(){
        this->a = p1;
        this->b = p2;
        useColors = false;
        this->amount = amount;
    }
};

Color patternAt(Pattern* pattern, Tuple point);

Color patternAt(Stripes* stripes, Tuple point){
    if(fmod(floor(point.x), 2) == 0){
        return stripes->a;
    }
    return stripes->b;
}

Color patternAt(Gradient* gradient, Tuple point){
    auto amount = point.x - floor(point.x);
    return interpolate(gradient->a, gradient->b, amount);
}

Color patternAt(Ring* ring, Tuple point){
    if(fmod(floor(sqrt((point.x * point.x) + (point.z * point.z))), 2) == 0){
        return ring->a;
    }
    return ring->b;
}

Color patternAt(Checker3D* checker, Tuple point){
    if(fmod(floor(point.x) + floor(point.y) + floor(point.z), 2) == 0){
        return checker->a;
    }
    return checker->b;
}

Color patternAt(NestedChecker3D* checker, Tuple point){
    if(fmod(floor(point.x) + floor(point.y) + floor(point.z), 2) == 0){
        return patternAt(&checker->a, point);
    }
    return patternAt(&checker->b, point);
}

Color patternAt(Blend* blend, Tuple point){
    return interpolate(patternAt(&blend->a, point), patternAt(&blend->b, point), blend->amount);
}

Color patternAt(Pattern* pattern, Tuple point){
    auto pStripes = dynamic_cast<Stripes*>(pattern);
    if(pStripes != nullptr){
        return patternAt(pStripes, point);
    }

    auto pGradient = dynamic_cast<Gradient*>(pattern);
    if(pGradient != nullptr){
        return patternAt(pGradient, point);
    }

    auto pRing = dynamic_cast<Ring*>(pattern);
    if(pRing != nullptr){
        return patternAt(pRing, point);
    }

    auto pChecker = dynamic_cast<Checker3D*>(pattern);
    if(pChecker != nullptr){
        return patternAt(pChecker, point);
    }

    auto pNestedChecker = dynamic_cast<NestedChecker3D*>(pattern);
    if(pNestedChecker != nullptr){
        return patternAt(pNestedChecker, point);
    }

    auto pBlend = dynamic_cast<Blend*>(pattern);
    if(pBlend != nullptr){
        return patternAt(pBlend, point);
    }

    return BLACK;
}

//TODO: Testar se as Patterns estão funcionando com Groups
Color patternAtObject(Pattern* pattern, Object *object, Tuple worldPoint){
    auto objectPoint = worldToObject(object, worldPoint);
    auto patternPoint = inverse(pattern->transform) * objectPoint;
    return patternAt(pattern, patternPoint);
}

#endif //RAYTRACERCHALLENGE_PATTERN_H
