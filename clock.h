//
// Created by Henrique on 06/06/2019.
//

#ifndef RAYTRACERCHALLENGE_CLOCK_H
#define RAYTRACERCHALLENGE_CLOCK_H

#include "transformations.h"
#include "canvas.h"

void drawPoint(Canvas* canvas, Tuple* point){
    auto x = (unsigned int) (canvas->width*0.5 + point->x);
    auto y = (unsigned int) (canvas->width*0.5 + point->z);

    writePixel(canvas, x, y, Color{1, 1, 1});
}

void drawClock(){
    auto c = canvas(250, 250);
    Tuple hourPoints[12];

    hourPoints[0] = point(0, 0, 100);
    drawPoint(&c, &hourPoints[0]);

    for(int i = 1; i < 12; i++){
        hourPoints[i] = hourPoints[i - 1] * rotationY(radians(360.f / 12));
        drawPoint(&c, &hourPoints[i]);
    }

    canvasToPNG(&c, "clock.png");
}

#endif //RAYTRACERCHALLENGE_CLOCK_H
