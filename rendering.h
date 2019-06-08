//
// Created by Henrique on 07/06/2019.
//

#ifndef RAYTRACERCHALLENGE_RENDERING_H
#define RAYTRACERCHALLENGE_RENDERING_H

#include "matrix.h"
#include "color.h"
#include "canvas.h"
#include <vector>
#include "ray.h"
#include "transformations.h"

struct Material{
    Color color{};
    float ambient;
    float diffuse;
    float specular;
    float shininess;

    Material(){
        color = Color{1, 1, 1};
        ambient = 0.1f;
        diffuse = 0.9f;
        specular = 0.9f;
        shininess = 200;
    }
};

struct Sphere {
    unsigned int id;
    Tuple origin{};
    float radius;
    Matrix4x4 transform;
    Material material;

    Sphere(Tuple o = point(0, 0, 0)){
        static unsigned int cont = 0;
        id = cont++;
        origin = o;
        radius = 1;
        transform = Matrix4x4();
        material = Material();
    }
};

struct Intersection{
    Sphere sphere;
    float time;
};

struct Computations{
    Intersection intersection;
    Tuple point;
    Tuple eyev;
    Tuple normalv;
    bool inside;
    Tuple overPoint;
};

struct Camera{
    unsigned int hSize;
    unsigned int vSize;
    float fieldOfView;
    float pixelSize;
    float halfWidth;
    float halfHeight;
    Matrix4x4 transform;

    Camera(unsigned int hSize, unsigned int vSize, float fieldOfView){
        this->hSize = hSize;
        this->vSize = vSize;
        this->fieldOfView = fieldOfView;
        this->transform = Matrix4x4();

        auto halfView = tanf(fieldOfView/2);
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
    std::vector<Sphere> spheres  = std::vector<Sphere>();;

    World() = default;
};

std::vector<Intersection> intersect(Ray ray, Sphere sphere){
    auto tRay = transform(ray, inverse(sphere.transform));

    auto sphere_to_ray = tRay.origin - sphere.origin;

    auto a = dot(tRay.direction, tRay.direction);
    auto b = 2 * dot(tRay.direction, sphere_to_ray);
    auto c = dot(sphere_to_ray, sphere_to_ray) - 1;

    auto discriminant = (b * b) - 4 * a * c;

    if(discriminant < 0){
        return std::vector<Intersection>();
    }

    auto raiz = sqrtf(discriminant);
    auto divisor = (2 * a);

    auto t1 = (-b - raiz) / divisor;
    auto t2 = (-b + raiz) / divisor;

    auto lista = std::vector<Intersection>();

    lista.push_back(Intersection{sphere, t1});
    lista.push_back(Intersection{sphere, t2});

    return lista;
}

Intersection hit(const std::vector<Intersection>& intersections){
    auto hit = Intersection{Sphere(), -1};
    for(auto intersection: intersections){
        if(intersection.time >= 0 && (intersection.time < hit.time || hit.time == -1)){
            hit = intersection;
        }
    }
    return hit;
}

Tuple normalAt(Sphere sphere, Tuple p){
    auto objectPoint = inverse(sphere.transform) * p;
    auto objectNormal = objectPoint - sphere.origin;
    auto worldNormal = transpose(inverse(sphere.transform)) * objectNormal;
    worldNormal.w = 0;
    return normalize(worldNormal);
}

Color lighting(Material material, Light light, Tuple point, Tuple eyev, Tuple normalv, bool inShadow){
    //Combine the surface color with the light's color/intensity
    auto effectiveColor = material.color * light.intensity;

    //Compute the ambient contribution
    auto ambient = effectiveColor * material.ambient;

    if(inShadow){
        return ambient;
    }

    //Find the direction to the light source
    auto lightv = normalize(light.position - point);

    //lightDotNormal represents the cosine of the angle between the light vector and the normal vector.
    //A negative number means the light is on the other side of the surface.
    auto lightDotNormal = dot(lightv, normalv);

    auto diffuse = Color{0, 0, 0};
    auto specular = Color{0, 0, 0};

    if(lightDotNormal >= 0){
        //Compute the diffuse contribuition
        diffuse = effectiveColor * material.diffuse * lightDotNormal;

        //reflectDotEye represents the cosine of the angle between the reflection vector and the eye vector.
        //A negative number means the light reflects away from the eye.
        auto reflectv = reflect(-lightv, normalv);
        auto reflectDotEye = dot(reflectv, eyev);

        if(reflectDotEye > 0){
            //Compute the specular contribuition
            auto factor = pow(reflectDotEye, material.shininess);
            specular = light.intensity * material.specular * factor;
        }
    }

    return ambient + diffuse + specular;
}

World defaultWorld(){
    auto world = World();

    auto defaulLight = pointLight(Tuple{-10, 10, -10}, Color{1, 1, 1});
    world.lightSources.push_back(defaulLight);

    auto s1 = Sphere();
    s1.material.color = Color{0.8, 1.0, 0.6};
    s1.material.diffuse = 0.7;
    s1.material.specular = 0.2;
    world.spheres.push_back(s1);

    auto s2 = Sphere();
    s2.transform = scaling(0.5, 0.5, 0.5);
    world.spheres.push_back(s2);

    return world;
}

std::vector<Intersection> intersectWorld(Ray ray, World world){
    auto lista = std::vector<Intersection>();
    for(const Sphere& sphere: world.spheres){
        auto intersecoes = intersect(ray, sphere);
        for(auto intersecao: intersecoes){
            unsigned int pos = 0;
            for(; pos < lista.size(); pos++){
                if(lista[pos].time > intersecao.time){
                    break;
                }
            }
            lista.insert(lista.begin() + pos, intersecao);
        }
    }
    return lista;
}

Computations prepareComputations(Ray ray, Intersection intersection){
    auto comps = Computations{};
    comps.intersection = intersection;
    comps.point = position(ray, intersection.time);
    comps.eyev = - ray.direction;
    comps.normalv = normalAt(intersection.sphere, comps.point);
    comps.inside = dot(comps.normalv, comps.eyev) < 0;

    if(comps.inside){
        comps.normalv = - comps.normalv;
    }

    comps.overPoint = comps.point + (comps.normalv * EPSILON);

    return comps;
}

bool isShadowed(World world, Tuple point){
    for(auto light: world.lightSources){
        auto v = light.position - point;
        auto distance = length(point);
        auto direction = normalize(v);

        auto r = Ray{ point, direction };
        auto intersections = intersectWorld(r, world);
        auto h = hit(intersections);
        if(h.time >= 0 && h.time < distance){
            return true;
        }
    }
    return false;
}

Color shadeHit(World world, Computations comps, bool shadows){
    auto material = comps.intersection.sphere.material;

    auto color = Color();
    for(auto light: world.lightSources){
        color = color + lighting(material, light, comps.overPoint, comps.eyev, comps.normalv, shadows && isShadowed(world, comps.overPoint));
    }
    return color;
}

Color colorAt(Ray ray, const World& world, bool shadows){
    auto intersections = intersectWorld(ray, world);
    auto h = hit(intersections);
    auto color = Color{0, 0, 0};

    if(h.time >= 0){
        auto comps = prepareComputations(ray, h);
        color = shadeHit(world, comps, shadows);
    }

    return color;
}

Matrix4x4 viewTransform(Tuple from, Tuple to, Tuple up){
    auto forward = normalize(to - from);
    auto left = cross(forward, normalize(up));
    auto trueUp = cross(left, forward);
    auto orientation = Matrix4x4(left.x,     left.y,     left.z,     0,
                                 trueUp.x,   trueUp.y,   trueUp.z,   0,
                                 -forward.x, -forward.y, -forward.z, 0,
                                 0,          0,          0,          1);
    return orientation * translation(-from.x, -from.y, -from.z);
}

Ray rayForPixel(Camera camera, unsigned int x, unsigned int y){
    auto xOffset = (x + 0.5f) * camera.pixelSize;
    auto yOffset = (y + 0.5f) * camera.pixelSize;

    auto worldX = camera.halfWidth - xOffset;
    auto worldY = camera.halfHeight - yOffset;

    auto pixel = inverse(camera.transform) * point(worldX, worldY, -1);
    auto origin = inverse(camera.transform) * point(0, 0, 0);
    auto direction = normalize(pixel - origin);

    return Ray{origin, direction};
}

Canvas render(Camera camera, const World& world, bool shadows){
    auto image = canvas(camera.hSize, camera.vSize);

    for(unsigned int y = 0; y < camera.vSize; y++){
        for(unsigned int x = 0; x < camera.hSize; x++){
            auto ray = rayForPixel(camera, x, y);
            auto color = colorAt(ray, world, shadows);
            writePixel(&image, x, y, color);
        }
    }

    return image;
}

void drawSphereRaycast(unsigned int canvasPixels){

    auto c = canvas(canvasPixels, canvasPixels);
    auto s = Sphere();

    auto rayOrigin = point(0, 0, -5);
    float wallSize = 7;
    float wallZ = 10;
    float pixelSize = wallSize / canvasPixels;
    float half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto position = point(worldX, worldY, wallZ);
            auto ray = Ray{rayOrigin, normalize(position - rayOrigin)};

            if(hit(intersect(ray, s)).time >= 0){
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
    float wallSize = 7;
    float wallZ = 10;
    float pixelSize = wallSize / canvasPixels;
    float half = wallSize / 2;

    for(unsigned int y = 0; y < c.height; y++){

        auto worldY = half - pixelSize * y;

        for(unsigned int x = 0; x < c.width; x++){

            auto worldX = - half + pixelSize * x;
            auto ray = Ray{rayOrigin, normalize(point(worldX, worldY, wallZ) - rayOrigin)};
            ray.direction = normalize(ray.direction);

            auto closest = hit(intersect(ray, s));

            if(closest.time >= 0){
                auto point = position(ray, closest.time);
                auto normal = normalAt(s, point);
                auto eye = - ray.direction;

                auto color = lighting(s.material, light, point, eye, normal, false);
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
    world.spheres.push_back(floor);
    world.spheres.push_back(leftWall);
    world.spheres.push_back(rightWall);
    world.spheres.push_back(middle);
    world.spheres.push_back(right);
    world.spheres.push_back(left);

    auto camera = Camera(width, height, PI/3);
    camera.transform = viewTransform(point(0, 1.5f, -5), point(0, 1, 0), vector(0, 1, 0));

    auto canvas = render(camera, world, shadows);
    canvasToPNG(&canvas, "world.png");
}

#endif //RAYTRACERCHALLENGE_RENDERING_H
