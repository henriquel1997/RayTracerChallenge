#include <chrono>
#include "render_functions.h"

int main() {
    using namespace std::chrono;
    auto startTime = high_resolution_clock::now();
    //drawSphereRaycast(100);
    //drawSpherePhong(100);
    //drawWorldScene(100, 50, true);
    //drawPlaneScene(100, 50, true);
    //drawPatternScene(100, 50, true);
    //drawReflectionScene(1000, 500, true);
    //drawGlass(1000, 1000, false);
    //drawGlassAndCube(1000, 100, false);
    drawCylinder(1000, 1000, true);
    auto endTime = high_resolution_clock::now();
    auto totalTime = duration_cast<milliseconds>(endTime - startTime).count();
    printf("Time: %lli(ms)\n", totalTime);
}