//
// Created by Henrique on 08/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RENDER_FUNCTIONS_H
#define RAYTRACERCHALLENGE_RENDER_FUNCTIONS_H

#include "rendering.h"

void drawSphereRaycast(unsigned int canvasPixels){

    auto c = canvas(canvasPixels, canvasPixels);
    auto s = Sphere();

    auto rayOrigin = point(0, 0, -5);
    double wallSize = 7;
    double wallZ = 10;
    double pixelSize = wallSize / canvasPixels;
    double half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto position = point(worldX, worldY, wallZ);
            auto ray = Ray{rayOrigin, normalize(position - rayOrigin)};

            if(hit(intersect(ray, &s)).time >= 0){
                writePixel(&c, x, y, Color{1, 0, 0});
            }
        }
    }

    canvasToPNG(&c, "sphere.png");
}

void drawSpherePhong(unsigned int canvasPixels){

    auto c = canvas(canvasPixels, canvasPixels);
    auto s = Sphere();
    s.material.color = Color{ 1, 0.2, 1 };

    auto light = pointLight(point(-10, 10, -10), Color{1, 1, 1});

    auto rayOrigin = point(0, 0, -5);
    double wallSize = 7;
    double wallZ = 10;
    double pixelSize = wallSize / canvasPixels;
    double half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto ray = Ray{rayOrigin, normalize(point(worldX, worldY, wallZ) - rayOrigin)};
            ray.direction = normalize(ray.direction);

            auto closest = hit(intersect(ray, &s));

            if(closest.time >= 0){
                auto point = position(ray, closest.time);
                auto normal = normalAt(&s, point);
                auto eye = - ray.direction;

                auto color = lighting(s.material, &s, light, point, eye, normal, false);
                writePixel(&c, x, y, color);
            }
        }
    }

    canvasToPNG(&c, "phong.png");
}

void drawWorldScene(unsigned int width, unsigned int height, bool shadows){
    auto floor = Sphere();
    floor.transform = scaling(10, 0.01, 10);
    floor.material.color = Color{ 1, 0.9f, 0.9f };
    floor.material.specular = 0;

    auto leftWall = Sphere();
    leftWall.transform = translation(0, 0, 5) * rotationY(-PI/4) * rotationX(PI/2) * scaling(10, 0.01, 10);
    leftWall.material = floor.material;

    auto rightWall = Sphere();
    rightWall.transform = translation(0, 0, 5) * rotationY(PI/4) * rotationX(PI/2) * scaling(10, 0.01, 10);
    rightWall.material = floor.material;

    auto middle = Sphere();
    middle.transform = translation(-.5f, 1.f, .5f);
    middle.material.color = Color{ 0.1f, 1, .5f };
    middle.material.diffuse = 0.7f;
    middle.material.specular = 0.3f;

    auto right = Sphere();
    right.transform = translation(1.5f, .5f, -.5f) * scaling(.5f, .5f, .5f);
    right.material.color = Color{ 0.5f, 1, .1f };
    right.material.diffuse = 0.7f;
    right.material.specular = 0.3f;

    auto left = Sphere();
    left.transform = translation(-1.5f, .33f, -.75f) * scaling(.33f, .33f, .33f);
    left.material.color = Color{ 1.f, .8f, .1f };
    left.material.diffuse = 0.7f;
    left.material.specular = 0.3f;

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&floor);
    world.objects.push_back(&leftWall);
    world.objects.push_back(&rightWall);
    world.objects.push_back(&middle);
    world.objects.push_back(&right);
    world.objects.push_back(&left);

    auto camera = Camera(width, height, PI/3);
    camera.transform = viewTransform(point(0, 1.5f, -5), point(0, 1, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "world.png");
}

void drawPlaneScene(unsigned int width, unsigned int height, bool shadows){

    auto plane = Plane();
    plane.material.color = Color{ 0.82f, 0.82f, 0.82f };
    plane.material.diffuse = 0.7f;
    plane.material.specular = 0.3f;

    auto middle = Sphere();
    middle.transform = translation(-.5f, 1.f, .5f);
    middle.material.color = Color{ 0.1f, 1, .5f };
    middle.material.diffuse = 0.7f;
    middle.material.specular = 0.3f;

    auto right = Sphere();
    right.transform = translation(1.5f, .5f, -.5f) * scaling(.5f, .5f, .5f);
    right.material.color = Color{ 0.5f, 1, .1f };
    right.material.diffuse = 0.7f;
    right.material.specular = 0.3f;

    auto left = Sphere();
    left.transform = translation(-1.5f, .33f, -.75f) * scaling(.33f, .33f, .33f);
    left.material.color = Color{ 1.f, .8f, .1f };
    left.material.diffuse = 0.7f;
    left.material.specular = 0.3f;

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&plane);
    world.objects.push_back(&middle);
    world.objects.push_back(&right);
    world.objects.push_back(&left);

    auto camera = Camera(width, height, PI/3);
    camera.transform = viewTransform(point(0, 1.5f, -5), point(0, 1, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "plane.png");
}

void drawPatternScene(unsigned int width, unsigned int height, bool shadows){
    auto c1 = GREEN;
    auto c2 = BLUE;
    auto c3 = RED;
    auto c4 = WHITE;
    auto c5 = BLACK;
    auto stripes = Stripes(c1, c2);
    auto ring = Ring(c3, c4);
    auto gradient = Gradient(c3, c2);
    auto checker = Checker3D(c5, c4);
    checker.transform = scaling(0.5, 0.5, 0.5);

    auto plane = Plane();
    plane.material.color = Color{ 0.82f, 0.82f, 0.82f };
    plane.material.diffuse = 0.7f;
    plane.material.specular = 0.3f;
    addPattern(&plane.material, &ring);

    auto middle = Sphere();
    middle.transform = translation(-.5f, 1.f, .5f);
    middle.material.color = Color{ 0.1f, 1, .5f };
    middle.material.diffuse = 0.7f;
    middle.material.specular = 0.3f;
    addPattern(&middle.material, &checker);

    auto right = Sphere();
    right.transform = translation(1.5f, .5f, -.5f) * scaling(.5f, .5f, .5f);
    right.material.color = Color{ 0.5f, 1, .1f };
    right.material.diffuse = 0.7f;
    right.material.specular = 0.3f;
    addPattern(&right.material, &stripes);

    auto left = Sphere();
    left.transform = translation(-1.5f, .33f, -.75f) * scaling(.33f, .33f, .33f);
    left.material.color = Color{ 1.f, .8f, .1f };
    left.material.diffuse = 0.7f;
    left.material.specular = 0.3f;
    addPattern(&left.material, &gradient);

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&plane);
    world.objects.push_back(&middle);
    world.objects.push_back(&right);
    world.objects.push_back(&left);

    auto camera = Camera(width, height, PI/3);
    camera.transform = viewTransform(point(0, 1.5f, -5), point(0, 1, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "pattern.png");
}

void drawReflectionScene(unsigned int width, unsigned int height, bool shadows){
    auto c1 = GREEN;
    auto c2 = BLUE;
    auto c3 = RED;
    auto c4 = WHITE;
    auto c5 = BLACK;

    auto stripes = Stripes(c1, c2);
    auto ring = Ring(c3, c4);
    auto gradient = Gradient(c3, c2);
    auto checker = Checker3D(c5, c4);
    checker.transform = scaling(0.5, 0.5, 0.5);

    auto plane = Plane();
    plane.material.color = Color{ 0.82f, 0.82f, 0.82f };
    plane.material.diffuse = 0.7f;
    plane.material.specular = 0.3f;
    plane.material.reflective = 1.f;
    addPattern(&plane.material, &ring);

    auto middle = Sphere();
    middle.transform = translation(-.5f, 1.f, .5f);
    middle.material.color = Color{ 0.1f, 1, .5f };
    middle.material.diffuse = 0.7f;
    middle.material.specular = 0.3f;
    middle.material.reflective = 1.f;
    addPattern(&middle.material, &checker);

    auto right = Sphere();
    right.transform = translation(1.5f, .5f, -.5f) * scaling(.5f, .5f, .5f);
    right.material.color = Color{ 0.5f, 1, .1f };
    right.material.diffuse = 0.7f;
    right.material.specular = 0.3f;
    right.material.reflective = 1.f;
    addPattern(&right.material, &stripes);

    auto left = Sphere();
    left.transform = translation(-1.5f, .33f, -.75f) * scaling(.33f, .33f, .33f);
    left.material.color = Color{ 1.f, .8f, .1f };
    left.material.diffuse = 0.7f;
    left.material.specular = 0.3f;
    left.material.reflective = 1.f;
    addPattern(&left.material, &gradient);

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&plane);
    world.objects.push_back(&middle);
    world.objects.push_back(&right);
    world.objects.push_back(&left);

    auto camera = Camera(width, height, PI/3);
    camera.transform = viewTransform(point(0, 1.5f, -5), point(0, 1, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "reflection.png");
}

Sphere glassSphere(){
    auto sphere = Sphere();
    sphere.material.transparency = 1.0;
    sphere.material.reflective = 1.0;
    sphere.material.refractiveIndex = 1.52;
    sphere.material.color = BLACK;
    return sphere;
}

Sphere airSphere(){
    auto sphere = Sphere();
    sphere.material.transparency = 1.0;
    sphere.material.reflective = 1.0;
    sphere.material.refractiveIndex = 1.00029;
    sphere.material.color = BLACK;
    return sphere;
}

void drawGlass(unsigned int width, unsigned int height, bool shadows){
    auto c1 = GREEN;
    auto c2 = BLUE;
    auto c3 = RED;
    auto c4 = WHITE;
    auto c5 = BLACK;

    auto stripes = Stripes(c1, c2);
    auto ring = Ring(c3, c4);
    auto gradient = Gradient(c3, c2);
    auto checker = Checker3D(c5, c4);
    checker.transform = scaling(0.1, 0.1, 0.1);

    auto plane = Plane();
    plane.material.color = WHITE;
    plane.material.diffuse = 0.7f;
    plane.material.specular = 0.3f;
    plane.transform = rotationX(PI/2);
    addPattern(&plane.material, &checker);

    auto glass = glassSphere();

    auto air = airSphere();
    air.transform = scaling(0.5, 0.5, 0.5);

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&plane);
    world.objects.push_back(&glass);
    world.objects.push_back(&air);

    auto camera = Camera(width, height, PI/6);
    camera.transform = viewTransform(point(0, 2, -3), point(0, 0, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "glass.png");
}

void drawGlassAndCube(unsigned int width, unsigned int height, bool shadows){
    auto c1 = GREEN;
    auto c2 = BLUE;
    auto c3 = RED;
    auto c4 = WHITE;
    auto c5 = BLACK;

    auto stripes = Stripes(c1, c2);
    auto ring = Ring(c3, c4);
    auto gradient = Gradient(c3, c2);
    auto checker = Checker3D(c5, c4);
    checker.transform = scaling(0.1, 0.1, 0.1);

    auto box = Cube();
    box.material.diffuse = 0.7f;
    box.material.specular = 0.3f;
    box.transform = scaling(10, 10, 10);
    addPattern(&box.material, &checker);

    auto glass = glassSphere();

    auto cube = Cube();
    cube.material.color = RED;
    cube.material.diffuse = 0.7f;
    cube.material.specular = 0.3f;
    cube.transform = scaling(0.1, 0.1, 0.1) * rotationY(PI/4);

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&box);
    world.objects.push_back(&glass);
    world.objects.push_back(&cube);

    auto camera = Camera(width, height, PI/6);
    camera.transform = viewTransform(point(0, 2, -3), point(0, 0, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "glass and cube.png");
}

void drawCylinder(unsigned int width, unsigned int height, bool shadows){
    auto c1 = GREEN;
    auto c2 = BLUE;
    auto c3 = RED;
    auto c4 = WHITE;
    auto c5 = BLACK;

    auto stripes = Stripes(c1, c2);
    auto ring = Ring(c3, c4);
    auto gradient = Gradient(c3, c2);
    auto checker = Checker3D(c5, c4);
    checker.transform = scaling(0.1, 0.1, 0.1);

    auto box = Cube();
    box.material.diffuse = 0.7f;
    box.material.specular = 0.3f;
    box.transform = scaling(10, 10, 10);
    addPattern(&box.material, &checker);

    auto cube = Cylinder();
    cube.minimumY = -1;
    cube. maximumY = 1;
    cube.closed = true;
    cube.material.color = RED;
    cube.material.diffuse = 0.7f;
    cube.material.specular = 0.3f;
    cube.transform = scaling(0.5, 0.5, 0.5) * rotationX(-PI/4) * rotationZ(-PI/4);

    auto world = World();
    world.lightSources.push_back(pointLight(point(-10, 10, -10), Color{ 1, 1, 1 }));
    world.objects.push_back(&box);
    world.objects.push_back(&cube);

    auto camera = Camera(width, height, PI/6);
    camera.transform = viewTransform(point(0, 1, -3), point(0, 0, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "cylinder.png");
}

#endif //RAYTRACERCHALLENGE_RENDER_FUNCTIONS_H
