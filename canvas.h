//
// Created by Henrique on 05/06/2019.
//

#ifndef RAYTRACERCHALLENGE_CANVAS_H
#define RAYTRACERCHALLENGE_CANVAS_H

#include <cstdlib>
#include "color.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

struct Canvas{
    unsigned int width;
    unsigned int height;
    Color* pixels;
};

Canvas canvas(unsigned int width, unsigned int height){
    auto pixels = (Color*) malloc(width * height * sizeof(Color));

    for(unsigned int i = 0; i < width * height; i++){
        pixels[i] = Color{ 0, 0, 0 };
    }

    return Canvas { width, height, pixels };
}

unsigned int getCanvasPos(Canvas* canvas, unsigned int x, unsigned int y){
    return x + (y * canvas->width);
}

Color* pixelAt(Canvas* canvas, unsigned int x, unsigned int y){
    if(x < canvas->width && y < canvas->height && canvas->pixels != nullptr){
        return canvas->pixels + getCanvasPos(canvas, x, y);
    }
    return nullptr;
}

void writePixel(Canvas* canvas, unsigned int x, unsigned int y, Color color){
    if(x < canvas->width && y < canvas->height && canvas->pixels != nullptr){
        auto pos = getCanvasPos(canvas, x, y);
        canvas->pixels[pos] = color;
    }
}

void canvasToPNG(Canvas* canvas, const char* pathAndName){

    auto canvasRGB  = new unsigned char [3 * canvas->width * canvas->height];

    for(unsigned int i = 0; i < canvas->width * canvas->height; i++){
        auto pixel = canvas->pixels[i];
        pixel.red = clamp(pixel.red, 0.f, 1.f);
        pixel.green = clamp(pixel.green, 0.f, 1.f);
        pixel.blue = clamp(pixel.blue, 0.f, 1.f);

        canvasRGB[(3*i) + 0] = (unsigned char) (255 * pixel.red);
        canvasRGB[(3*i) + 1] = (unsigned char) (255 * pixel.green);
        canvasRGB[(3*i) + 2] = (unsigned char) (255 * pixel.blue);
    }

    stbi_write_png(pathAndName, canvas->width, canvas->height, 3, canvasRGB, 3 * canvas->width);
}

#endif //RAYTRACERCHALLENGE_CANVAS_H
