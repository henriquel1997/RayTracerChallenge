#include <cstdio>
#include "projectile.h"
#include "matrix.h"

int main() {
    //executeTick();

    auto matrix = Matrix4x4();

    matrix[0][0] = 10;

    printf("[0][0]: %f", matrix[0][0]);

    return 0;
}