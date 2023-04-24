#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

    // TODO Project 3-2: Part 4
    // compute position and direction of ray from the input sensor sample coordinate.
    // Note: use rndR and rndTheta to uniformly sample a unit disk.
    double sensor_width = 2.0 * tan(hFov * PI * 0.5 / 180);
    double sensor_height = 2.0 * tan(vFov * PI * 0.5 / 180);

    // Calculate position of sensor sample on sensor plane
    double sensor_x = sensor_width * (x - 0.5);
    double sensor_y = sensor_height * (y - 0.5);

    Vector3D sensor_point(sensor_x, sensor_y, -1.0);
    // transform the sensor point to the camera coordinate system
    sensor_point = (c2w * sensor_point).unit();

    // Uniformly sample the disk representing the thin lens
    Vector3D pLens = lensRadius * sqrt(rndR) * Vector3D(cos(rndTheta), sin(rndTheta), 0);

    // Calculate the position and direction of the ray from the lens to the focus plane
    Vector3D pFocus = sensor_point * focalDistance;
    Vector3D direction = (pFocus - pLens).unit();
    Vector3D origin = c2w * pLens + pos;

    Ray ray = Ray(origin, c2w * direction);
    ray.max_t = fClip;
    ray.min_t = nClip;
    return ray;
}


} // namespace CGL
