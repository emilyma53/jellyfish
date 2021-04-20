#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass& pm) {
	// TODO (Part 3): Handle collisions with spheres.
	if ((pm.position - origin).norm() <= radius) {
		Vector3D tan = (pm.position - origin).unit() * radius + origin;
		Vector3D correct = tan - pm.last_position;
		pm.position = (1 - friction) * correct + pm.last_position;
	}
}

void Sphere::render(GLShader& shader) {
	// We decrease the radius here so flat triangles don't behave strangely
	// and intersect with the sphere when rendered
	m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
