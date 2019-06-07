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

struct Object{
    unsigned int id;
    Tuple origin{};

    explicit Object(Tuple o = point(0, 0, 0)){
        static unsigned int cont = 0;
        id = cont++;
        origin = o;
    }

    virtual ~Object() = default;
};

struct Sphere: Object {
    float radius;
    Matrix4x4 transform;
    Material material;

    Sphere(){
        radius = 1;
        transform = Matrix4x4();
        material = Material();
    }
};

struct Intersection{
    Object* object;
    float time;
};

struct Computations{
    Intersection intersection;
    Tuple point;
    Tuple eyev;
    Tuple normalv;
    bool inside;
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
    std::vector<Object> objects  = std::vector<Object>();;

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

    lista.push_back(Intersection{&sphere, t1});
    lista.push_back(Intersection{&sphere, t2});

    return lista;
}

std::vector<Intersection> intersect(Ray ray, Object object){
    auto pSphere = dynamic_cast<Sphere*>(&object);
    if(pSphere != nullptr){
        return intersect(ray, *pSphere);
    }
    return std::vector<Intersection>();
}

Intersection hit(const std::vector<Intersection>& intersections){
    auto hit = Intersection{nullptr, -1};
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

Tuple normalAt(Object object, Tuple p){
    auto pSphere = dynamic_cast<Sphere*>(&object);
    if(pSphere != nullptr){
        return normalAt(*pSphere, p);
    }
    return vector(0, 0, 0);
}

Color lighting(Material material, Light light, Tuple point, Tuple eyev, Tuple normalv){
    //Combine the surface color with the light's color/intensity
    auto effectiveColor = material.color * light.intensity;

    //Find the direction to the light source
    auto lightv = normalize(light.position - point);

    //Compute the ambient contribution
    auto ambient = effectiveColor * material.ambient;

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
    world.objects.push_back(s1);

    auto s2 = Sphere();
    s2.transform = scaling(0.5, 0.5, 0.5);
    world.objects.push_back(s2);

    return world;
}

std::vector<Intersection> intersectWorld(Ray ray, World world){
    auto lista = std::vector<Intersection>();
    for(const Object& object: world.objects){
        for(auto intersecao: intersect(ray, object)){
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
    comps.normalv = normalAt(*intersection.object, comps.point);
    comps.inside = dot(comps.normalv, comps.eyev) < 0;

    if(comps.inside){
        comps.normalv = - comps.normalv;
    }

    return comps;
}

Color shadeHit(World world, Computations comps){
    auto material = Material();
    auto pSphere = dynamic_cast<Sphere*>(comps.intersection.object);
    if(pSphere != nullptr){
        material = pSphere->material;
    }

    auto color = Color();
    for(auto light: world.lightSources){
        color = color + lighting(material, light, comps.point, comps.eyev, comps.normalv);
    }
    return color;
}

Color colorAt(Ray ray, const World& world){
    auto intersections = intersectWorld(ray, world);
    auto h = hit(intersections);
    auto color = Color{0, 0, 0};

    if(h.time >= 0){
        auto comps = prepareComputations(ray, h);
        color = shadeHit(world, comps);
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

Canvas render(Camera camera, const World& world){
    auto image = canvas(camera.hSize, camera.vSize);

    for(unsigned int y = 0; y < camera.vSize; y++){
        for(unsigned int x = 0; x < camera.hSize; x++){
            auto ray = rayForPixel(camera, x, y);
            auto color = colorAt(ray, world);
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

                auto color = lighting(s.material, light, point, eye, normal);
                writePixel(&c, x, y, color);
            }
        }
    }

    canvasToPNG(&c, "phong.png");
}

#endif //RAYTRACERCHALLENGE_RENDERING_H
