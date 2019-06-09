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
#include "pattern.h"
#include "structs.h"

#define MAX_DEPTH 5
Color colorAt(Ray ray, const World& world, bool shadows, unsigned int remaining = MAX_DEPTH);
Color reflectedColor(const World &world, Computations comps, bool shadows, unsigned int remaining);

void addPattern(Material* material, Pattern* pattern){
    material->hasPattern = true;
    material->pattern = pattern;
}

Intersection hit(const std::vector<Intersection>& intersections, bool ignoreNoShadows = false){
    auto hit = Intersection{nullptr, -1};
    for(const auto &intersection: intersections){
        //When ignoreNoShadows is true, the function will ignore objects that doesn't cast shadows
        if(ignoreNoShadows && !intersection.object->material.castShadows){
            continue;
        }

        if(intersection.time >= 0 && (intersection.time < hit.time || hit.time == -1)){
            hit = intersection;
        }
    }
    return hit;
}

Color lighting(Material material, Object* object, Light light, Tuple point, Tuple eyev, Tuple normalv, bool inShadow){

    Color color{};
    if(material.hasPattern){
        color = patternAtObject(material.pattern, object, point);
    }else{
        color = material.color;
    }

    //Combine the surface color with the light's color/intensity
    auto effectiveColor = color * light.intensity;

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

    auto diffuse = BLACK;
    auto specular = BLACK;

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

    auto color = WHITE;
    auto defaulLight = pointLight(Tuple{-10, 10, -10}, color);
    world.lightSources.push_back(defaulLight);

    auto s1 = Sphere();
    s1.material.color = Color{0.8, 1.0, 0.6};
    s1.material.diffuse = 0.7;
    s1.material.specular = 0.2;
    world.objects.push_back(&s1);

    auto s2 = Sphere();
    s2.transform = scaling(0.5, 0.5, 0.5);
    world.objects.push_back(&s2);

    return world;
}

std::vector<Intersection> intersectWorld(Ray ray, World world){
    auto lista = std::vector<Intersection>();
    for(Object* object: world.objects){
        auto intersecoes = intersect(ray, object);
        for(const auto &intersecao: intersecoes){
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

bool operator == (Intersection i1, Intersection i2){
    return i1.time == i2.time && ((i1.object == nullptr && i2.object == nullptr) || (i1.object->id == i2.object->id));
}

bool operator == (Object o1, Object o2){
    return o1.id == o2.id;
}

int contains(std::vector<Object> objects, Object* o){
    int cont = 0;
    for(const auto &item: objects){
        if(item == *o){
            return cont;
        }
        cont++;
    }
    return -1;
}

Computations prepareComputations(Ray ray, Intersection hit, std::vector<Intersection> intersections){
    auto comps = Computations{};
    comps.intersection = hit;
    comps.point = position(ray, hit.time);
    comps.eyev = - ray.direction;
    comps.normalv = normalAt(hit.object, comps.point);
    comps.inside = dot(comps.normalv, comps.eyev) < 0;

    if(comps.inside){
        comps.normalv = - comps.normalv;
    }

    auto smallNormal = (comps.normalv * EPSILON);
    comps.overPoint = comps.point + smallNormal;
    comps.underPoint = comps.point - smallNormal;

    //Refraction
    comps.reflectv = reflect(ray.direction, comps.normalv);

    auto containers = std::vector<Object>();
    for(auto intersection: intersections){
        if(intersection == hit){
            if(containers.empty()){
                comps.n1 = 1;
            }else{
                comps.n1 = containers[containers.size()-1].material.refractiveIndex;
            }
        }

        auto pos = contains(containers, intersection.object);
        if(pos >= 0){
            containers.erase(containers.begin() + pos);
        }else{
            containers.push_back(*intersection.object);
        }

        if(intersection == hit){
            if(containers.empty()){
                comps.n2 = 0;
            }else{
                comps.n2 = containers[containers.size()-1].material.refractiveIndex;
            }
        }
    }

    return comps;
}

bool isShadowed(World world, Tuple point){
    for(auto light: world.lightSources){
        auto v = light.position - point;
        auto distance = length(point);
        auto direction = normalize(v);

        auto r = Ray{ point, direction };
        auto intersections = intersectWorld(r, world);
        //Only want to check hits with objects that cast shadows
        auto h = hit(intersections, true);
        if(h.time >= 0 && h.time < distance){
            return true;
        }
    }
    return false;
}

Color refractedColor(const World &world, Computations comps, bool shadows, unsigned int remaining){
    auto transparency = comps.intersection.object->material.transparency;

    if(transparency == 0 || remaining <= 0){
        return BLACK;
    }

    auto nRatio = comps.n1 / comps.n2;
    auto cosI = dot(comps.eyev, comps.normalv);
    auto sin2T = nRatio*nRatio * (1 - (cosI*cosI));

    if(sin2T > 1){
        return BLACK;
    }

    auto cosT = sqrt(1 - sin2T);
    auto direction = (comps.normalv * ((nRatio * cosI) - cosT)) - (comps.eyev * nRatio);

    auto refractRay = Ray{ comps.underPoint, direction };

    return colorAt(refractRay, world, shadows, remaining -1) * transparency;
}

double schlick(Computations comps){
    auto cos = dot(comps.eyev, comps.normalv);

    if(comps.n1 > comps.n2){
        auto n = comps.n1 / comps.n2;
        auto sin2T = n*n * (1 - (cos*cos));

        if(sin2T > 1){
            return 1;
        }

        cos = sqrt(1 - sin2T);
    }

    auto sqrt_r0 = (comps.n1 - comps.n2) / (comps.n1 + comps.n2);
    auto r0 = sqrt_r0 * sqrt_r0;
    auto compCos = 1 - cos;
    auto compCos5 = compCos*compCos*compCos*compCos*compCos;
    return r0 + (1 - r0) * compCos5;
}

Color shadeHit(World world, Computations comps, bool shadows, unsigned int remaining){
    auto material = comps.intersection.object->material;

    auto surface = BLACK;
    for(auto light: world.lightSources){
        auto shadowed = shadows && isShadowed(world, comps.overPoint);
        surface = surface + lighting(material, comps.intersection.object, light, comps.overPoint, comps.eyev, comps.normalv, shadowed);
    }

    auto reflected = reflectedColor(world, comps, shadows, remaining);
    auto refracted = refractedColor(world, comps, shadows, remaining);

    if(material.reflective > 0 && material.transparency > 0){
        auto reflectance = schlick(comps);
        return surface + (reflected * reflectance) + (refracted * (1 - reflectance));
    }

    return surface + reflected + refracted;
}


Color colorAt(Ray ray, const World& world, bool shadows, unsigned int remaining){
    auto intersections = intersectWorld(ray, world);
    auto h = hit(intersections);
    auto color = BLACK;

    if(h.time >= 0){
        auto comps = prepareComputations(ray, h, intersections);
        color = shadeHit(world, comps, shadows, remaining);
    }

    return color;
}

Color reflectedColor(const World &world, Computations comps, bool shadows, unsigned int remaining){
    auto reflective = comps.intersection.object->material.reflective;
    if(remaining <= 0 || reflective <= 0){
        return BLACK;
    }

    auto reflectRay = Ray{ comps.overPoint, comps.reflectv };
    auto color = colorAt(reflectRay, world, shadows, remaining - 1);

    return color * reflective;
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

#endif //RAYTRACERCHALLENGE_RENDERING_H
