#include <cstdio>
#include "tuples.h"

int main() {
    auto p = point(4, -4, 3);
    printf("Point: ");
    printTuple(&p);
    printf("\n");

    auto v = vector(4, -4, 3);
    printf("Vector: ");
    printTuple(&v);
    printf("\n");

    printf("Iguais: %i", p == v);
    printf("\n");

    printf("Soma: ");
    auto soma = (p + v);
    printTuple(&soma);
    printf("\n");

    printf("Subtracao: ");
    auto sub = (p - v);
    printTuple(&sub);
    printf("\n");

    printf("Negacao: ");
    auto neg = (-p);
    printTuple(&neg);
    printf("\n");

    printf("Multiplicacao: ");
    auto mult = (3.5 * p);
    printTuple(&mult);
    printf("\n");

    printf("Divisao: ");
    auto div = (p / 2);
    printTuple(&div);
    printf("\n");

    printf("Magnitude: %f\n", length(&p));
    printf("Magnitude Quadrada: %f\n", lengthSquared(&p));

    printf("Normalizacao: ");
    auto norm = normalize(&p);
    printTuple(&norm);
    printf("\n");
    printf("Magnitude da Normalizacao: %f\n", length(&norm));

    printf("Dot Product: %f\n", dot(&p, &p));

    return 0;
}