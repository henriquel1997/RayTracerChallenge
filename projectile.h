//
// Created by Henrique on 05/06/2019.
//

#ifndef RAYTRACERCHALLENGE_PROJECTILE_H
#define RAYTRACERCHALLENGE_PROJECTILE_H

#include "tuples.h"
#include "canvas.h"

struct Projectile{
    Tuple position;
    Tuple velocity;
};

struct Environment{
    Tuple gravity;
    Tuple wind;
};

void tick(Environment* env, Projectile* proj){
    proj->position = proj->position + proj->velocity;
    proj->velocity = proj->velocity + env->gravity + env->wind;
}

void drawProjectile(Canvas* canvas, Projectile* projectile){
    auto x = (unsigned int) projectile->position.x;
    auto y = (unsigned int) (canvas->height - projectile->position.y);

    if(x >= 0 && x < canvas->width && y >= 0 && y < canvas->height){
        writePixel(canvas, x, y, Color{ 1, 0, 0 });
    }
}

void executeTick(){
    Projectile p = Projectile{point(0, 1, 0), normalize(vector(1, 1.8, 0)) * 11.25 };
    Environment e = Environment { vector(0, -0.1f, 0), vector(-0.01f, 0, 0) };

    auto c = canvas(900, 550);

    drawProjectile(&c, &p);
    while(p.position.y > 0){
        tick(&e, &p);
        drawProjectile(&c, &p);
    }

    canvasToPNG(&c, "projectile.png");
}

#endif //RAYTRACERCHALLENGE_PROJECTILE_H
