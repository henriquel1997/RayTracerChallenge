//
// Created by Henrique on 05/06/2019.
//

#ifndef RAYTRACERCHALLENGE_COLOR_H
#define RAYTRACERCHALLENGE_COLOR_H

#include "tuples.h"

struct Color{
    float red;
    float green;
    float blue;
};

Color add(Color* c1, Color* c2){
    return Color{ c1->red + c2->red, c1->green + c2->green, c1->blue + c2->blue };
}

Color operator + (Color c1, Color c2){
    return add(&c1, &c2);
}

Color subtract(Color* c1, Color* c2){
    return Color{ c1->red - c2->red, c1->green - c2->green, c1->blue - c2->blue };
}

Color operator - (Color c1, Color c2){
    return subtract(&c1, &c2);
}

Color multByScalar(Color* c, float scalar){
    return Color{ c->red * scalar, c->green * scalar, c->blue * scalar };
}

Color operator * (Color c, float scalar){
    return multByScalar(&c, scalar);
}

Color operator * (float scalar, Color c){
    return multByScalar(&c, scalar);
}

Color multByColor(Color* c1, Color* c2){
    return Color{ c1->red * c2->red, c1->green * c2->green, c1->blue * c2->blue };
}

Color operator * (Color c1, Color c2){
    return multByColor(&c1, &c2);
}

#endif //RAYTRACERCHALLENGE_COLOR_H
