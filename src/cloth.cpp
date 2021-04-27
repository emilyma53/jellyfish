#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
	int num_height_points, float thickness) {
	this->width = width;
	this->height = height;
	this->num_width_points = num_width_points;
	this->num_height_points = num_height_points;
	this->thickness = thickness;

	buildGrid();
	buildClothMesh();
}

Cloth::~Cloth() {
	point_masses.clear();
	springs.clear();

	if (clothMesh) {
		delete clothMesh;
	}
}

void Cloth::buildGrid() {
	// TODO (Part 1): Build a grid of masses and springs.

	// Add all point masses
//	double height_interval = this->height / this->num_height_points;
//	double width_interval = this->width / this->num_width_points;
//	for (int i = 0; i < this->num_height_points; i++) {
//		for (int j = 0; j < this->num_width_points; j++) {
//			Vector3D pos;
//			bool pin_point = false;
//			if (this->orientation == HORIZONTAL) {
//				pos = Vector3D(j * width_interval, 1.0f, i * height_interval);
//			}
//			else {
//				float z_pos = (rand() % 200) / 100000.0f - 0.001f;
//				pos = Vector3D(j * width_interval, i * height_interval, z_pos);
//			}
//			for (vector<int>& xy : this->pinned) {
//				if (xy[0] == i && xy[1] == j) {
//					pin_point = true;
//				}
//			}
//			this->point_masses.emplace_back(PointMass(pos, pin_point));
//		}
//	}
    num_width_points = 20;
    num_height_points = 7;
    this->point_masses = std::vector<PointMass>();
    this->springs = std::vector<Spring>();
    std::vector<double> T {20.0,20.0,20.0,20.0,20.0,20.0,20.0};
        std::vector<double> R {.0001, .5, 1.5, 2.5, 3.5, 4.5, 5.5};
        
        for (int i = 0; i < R.size(); i++) {
            for (int j = 0; j < T[i]; j++) {
                double r = R[i];
                double theta = double(j) * (2.0 * PI / T[i]);
                double x = r * cos(theta);
                double y = r * sin(theta);
                double z = -.1*(R[i] * R[i]);
//                double z = 0.0;
                Vector3D pos = Vector3D(x, y, z);
                this->point_masses.emplace_back(PointMass(pos, false));
            }
        }
    
    
        int counter = -1;
        for (int i = 0; i < R.size(); i++) {
            for (int j = 0; j < T[i]; j++) {
                counter++;
                if (j == 0) continue;
                PointMass* o = &this->point_masses[counter - 1];
                PointMass* p = &this->point_masses[counter];
                this->springs.emplace_back(Spring(o, p, STRUCTURAL));
                if (j + 1 == int(T[i])) {
                    o = &this->point_masses[counter - (19)];
                    this->springs.emplace_back(Spring(o, p, STRUCTURAL));
                }
            }
        }

    // add falange point masses
   
    Vector3D pos = point_masses[counter - 20].position;
    pos.z = pos.z - 4.0;
    this->point_masses.emplace_back(PointMass(pos, false));
    Vector3D pos2 = point_masses[counter - 19].position;
    pos2.z = pos2.z - 4.0;
    this->point_masses.emplace_back(PointMass(pos2, false));
    PointMass* o = &this->point_masses[counter - 20];
    PointMass* p = &this->point_masses[counter];
    this->springs.emplace_back(Spring(o, p, STRUCTURAL));
    o = &this->point_masses[counter - 19];
    p = &this->point_masses[counter + 1];
    this->springs.emplace_back(Spring(o, p, STRUCTURAL));
    

	// Add all springs
	// 1. Structural constraints exist between a point mass and the point mass to its left as well as the point mass above it.
	// 2. Shearing constraints exist between a point mass and the point mass to its diagonal upper left as well as the point mass to its diagonal upper right.
	// 3. Bending constraints exist between a point mass and the point mass two away to its left as well as the point mass two above it.
//	for (int i = 0; i < this->num_height_points; i++) {
//		for (int j = 0; j < this->num_width_points; j++) {
//			PointMass* p = &this->point_masses[i * this->num_width_points + j];
//			PointMass* o;
//			// structural constraints
//			if (j > 0) { // left
//				o = &this->point_masses[i * this->num_width_points + (j - 1)];
//				this->springs.emplace_back(Spring(o, p, STRUCTURAL));
//			}
//			if (i > 0) { // top
//				o = &this->point_masses[(i - 1) * this->num_width_points + j];
//				this->springs.emplace_back(Spring(o, p, STRUCTURAL));
//			}
//			// shearing constraints
//			if (i > 0 && j > 0) { // diagonal left
//				o = &this->point_masses[(i - 1) * this->num_width_points + (j - 1)];
//				this->springs.emplace_back(Spring(o, p, SHEARING));
//			}
//			if (i > 0 && j < (this->num_width_points - 1)) { // diagonal right
//				o = &this->point_masses[(i - 1) * this->num_width_points + (j + 1)];
//				this->springs.emplace_back(Spring(o, p, SHEARING));
//			}
//			// bending constraints
//			if (j > 1) { // left
//				o = &this->point_masses[i * this->num_width_points + (j - 2)];
//				this->springs.emplace_back(Spring(o, p, BENDING));
//			}
//			if (i > 1) { // top
//				o = &this->point_masses[(i - 2) * this->num_width_points + j];
//				this->springs.emplace_back(Spring(o, p, BENDING));
//			}
//		}
//	}
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters* cp,
	vector<Vector3D> external_accelerations,
	vector<CollisionObject*>* collision_objects) {
	double mass = width * height * cp->density / num_width_points / num_height_points;
	double delta_t = 1.0f / frames_per_sec / simulation_steps;

	// TODO (Part 2): Compute total force acting on each point mass.
	// sum external forces on each mass
	for (int i = 0; i < this->point_masses.size(); i++) {
		PointMass* pm = &this->point_masses[i];
		pm->forces = Vector3D();
		for (const Vector3D &a : external_accelerations) {
			pm->forces += mass * a;
		}
	}

	// apply spring forces to each mass
	for (int i = 0; i < this->springs.size(); i++) {
		Spring* s = &this->springs[i];
		if (s->spring_type == STRUCTURAL && cp->enable_structural_constraints) {
			Vector3D dp = s->pm_b->position - s->pm_a->position;
			Vector3D fs = cp->ks * (dp.norm() - s->rest_length) * dp.unit();
			s->pm_a->forces += fs;
			s->pm_b->forces -= fs;
		}
		else if (s->spring_type == SHEARING && cp->enable_shearing_constraints) {
			Vector3D dp = s->pm_b->position - s->pm_a->position;
			Vector3D fs = cp->ks * (dp.norm() - s->rest_length) * dp.unit();
			s->pm_a->forces += fs;
			s->pm_b->forces -= fs;
		}
		else if (s->spring_type == BENDING && cp->enable_bending_constraints) {
			Vector3D dp = s->pm_b->position - s->pm_a->position;
			Vector3D fs = 0.2 * cp->ks * (dp.norm() - s->rest_length) * dp.unit();
			s->pm_a->forces += fs;
			s->pm_b->forces -= fs;
		}
	}

	// TODO (Part 2): Use Verlet integration to compute new point mass positions
	for (int i = 0; i < this->point_masses.size(); i++) {
		if (!this->point_masses[i].pinned) {
			PointMass* pm = &this->point_masses[i];
			Vector3D dxdt = (1 - cp->damping / 100.) * (pm->position - pm->last_position) + pm->forces / mass * delta_t * delta_t;
			pm->last_position = pm->position;
			pm->position += dxdt;
		}
	}

	// TODO (Part 4): Handle self-collisions.
	build_spatial_map();
	for (int i = 0; i < this->point_masses.size(); i++) {
		self_collide(this->point_masses[i], simulation_steps);
	}

	// TODO (Part 3): Handle collisions with other primitives.
	for (int i = 0; i < this->point_masses.size(); i++) {
		for (CollisionObject* o : *collision_objects) {
			o->collide(this->point_masses[i]);
		}
	}

	// TODO (Part 2): Constrain the changes to be such that the spring does not change
	// in length more than 10% per timestep [Provot 1995].
	for (int i = 0; i < this->springs.size(); i++) {
		Spring* s = &this->springs[i];
		Vector3D dp = s->pm_b->position - s->pm_a->position;
		double length = dp.norm();
		if (length > s->rest_length * 1.1) {
			Vector3D maxspring = s->rest_length * 1.1 * dp.unit();
			if (s->pm_a->pinned) {
				s->pm_b->position = s->pm_a->position + maxspring;
			}
			else if (s->pm_b->pinned) {
				s->pm_a->position = s->pm_b->position - maxspring;
			}
			else {
				Vector3D midpt = (s->pm_a->position + s->pm_b->position) / 2;
				s->pm_a->position = midpt - maxspring / 2;
				s->pm_b->position = midpt + maxspring / 2;
			}
		}
	}
}

void Cloth::build_spatial_map() {
	for (const auto& entry : map) {
		delete(entry.second);
	}
	this->map.clear();
    // TODO (Part 4): Build a spatial map out of all of the point masses.
	for (PointMass &p : this->point_masses) {
	    float key = hash_position(p.position);
        if (this->map.find(key) == this->map.end()) {
            auto *bucket = new vector<PointMass *>();
            this->map[key] = bucket;
        }
        this->map[key]->push_back(&p);
    }
}

void Cloth::self_collide(PointMass& pm, double simulation_steps) {
	// TODO (Part 4): Handle self-collision for a given point mass.

	Vector3D avgcorrect = Vector3D();
	vector<PointMass*>* bucket = this->map[hash_position(pm.position)];
	int count = 0;
	for (PointMass* candidate : *bucket) {
		if (!(pm.position == candidate->position)) {
			Vector3D dist = pm.position - candidate->position;
			if (dist.norm() < 2 * this->thickness) {
				Vector3D correct = (2 * this->thickness - dist.norm()) * dist.unit();
				avgcorrect += correct;
				count += 1;
			}
		}
	}
	avgcorrect /= max(count, 1);
	pm.position += avgcorrect / simulation_steps;
}

float Cloth::hash_position(Vector3D pos) {
	// TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    float w = 3.0f * width / num_width_points;
    int x = floor(pos.x / w);
    float h = 3.0f * height / num_height_points;
    int y = floor(pos.y / h);
    float t = std::max(w, h);
    int z = floor(pos.z / t);
    return x << 10 ^ y << 5 ^ z;
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
	PointMass* pm = &point_masses[0];
	for (int i = 0; i < point_masses.size(); i++) {
		pm->position = pm->start_position;
		pm->last_position = pm->start_position;
		pm++;
	}
}

void Cloth::buildClothMesh() {
	if (point_masses.size() == 0) return;

	ClothMesh* clothMesh = new ClothMesh();
	vector<Triangle*> triangles;

	// Create vector of triangles
	for (int y = 0; y < num_height_points - 1; y++) {
		for (int x = 0; x < num_width_points; x++) {
			PointMass* pm = &point_masses[y * num_width_points + x];
			// Get neighboring point masses:
			/*                      *
			 * pm_A -------- pm_B   *
			 *             /        *
			 *  |         /   |     *
			 *  |        /    |     *
			 *  |       /     |     *
			 *  |      /      |     *
			 *  |     /       |     *
			 *  |    /        |     *
			 *      /               *
			 * pm_C -------- pm_D   *
			 *                      *
			 */

//			float u_min = x;
//			u_min /= num_width_points - 1;
//			float u_max = x + 1;
//			u_max /= num_width_points - 1;
            float u_min = x;
            u_min /= num_width_points;
            float u_max = x + 1;
            u_max /= num_width_points;
            float v_min = y;
            v_min /= num_height_points;
            float v_max = y + 1;
            v_max /= num_height_points;
            PointMass* pm_A = pm;
            PointMass* pm_B = pm + 1;
            PointMass* pm_C = pm + num_width_points;
            PointMass* pm_D = pm + num_width_points + 1;
            if (x + 1 == num_width_points) {
                PointMass* pm_A = pm;
                PointMass* pm_B = &point_masses[y * num_width_points];
                PointMass* pm_C = pm + num_width_points;
                PointMass* pm_D = pm_B + num_width_points;
                
            }

			Vector3D uv_A = Vector3D(u_min, v_min, 0);
			Vector3D uv_B = Vector3D(u_max, v_min, 0);
			Vector3D uv_C = Vector3D(u_min, v_max, 0);
			Vector3D uv_D = Vector3D(u_max, v_max, 0);


			// Both triangles defined by vertices in counter-clockwise orientation
			triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
				uv_A, uv_C, uv_B));
			triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
				uv_B, uv_C, uv_D));
		}
	}
//    for (int x = 0; x < num_width_points; x++) {
//        PointMass* pm = &point_masses[num_width_points + x];
//        PointMass* pm_C = pm;
//        PointMass* pm_B = &point_masses[num_width_points * num_height_points];
//        PointMass* pm_D = pm + 1;
//        float u_min = x;
//        u_min /= num_width_points;
//        float u_max = x + 1;
//        u_max /= num_width_points;
//        float v_min = 0.0;
//        v_min /= num_height_points;
//        float v_max = 1;
//        v_max /= num_height_points;
////        Vector3D uv_A = Vector3D(u_min, v_min, 0);
//        Vector3D uv_B = Vector3D(u_max, v_min, 0);
//        Vector3D uv_C = Vector3D(u_min, v_max, 0);
//        Vector3D uv_D = Vector3D(u_max, v_max, 0);
//        triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
//            uv_B, uv_C, uv_D));
//
//    }

	// For each triangle in row-order, create 3 edges and 3 internal halfedges
	for (int i = 0; i < triangles.size(); i++) {
		Triangle* t = triangles[i];

		// Allocate new halfedges on heap
		Halfedge* h1 = new Halfedge();
		Halfedge* h2 = new Halfedge();
		Halfedge* h3 = new Halfedge();

		// Allocate new edges on heap
		Edge* e1 = new Edge();
		Edge* e2 = new Edge();
		Edge* e3 = new Edge();

		// Assign a halfedge pointer to the triangle
		t->halfedge = h1;

		// Assign halfedge pointers to point masses
		t->pm1->halfedge = h1;
		t->pm2->halfedge = h2;
		t->pm3->halfedge = h3;

		// Update all halfedge pointers
		h1->edge = e1;
		h1->next = h2;
		h1->pm = t->pm1;
		h1->triangle = t;

		h2->edge = e2;
		h2->next = h3;
		h2->pm = t->pm2;
		h2->triangle = t;

		h3->edge = e3;
		h3->next = h1;
		h3->pm = t->pm3;
		h3->triangle = t;
	}

	// Go back through the cloth mesh and link triangles together using halfedge
	// twin pointers

	// Convenient variables for math
//	int num_height_tris = (num_height_points - 1) * 2;
//	int num_width_tris = (num_width_points - 1) * 2;
    int num_height_tris = (num_height_points - 1) * 2;
    int num_width_tris = (num_width_points) * 2;


	bool topLeft = true;
	for (int i = 0; i < triangles.size(); i++) {
		Triangle* t = triangles[i];

		if (topLeft) {
			// Get left triangle, if it exists
			if (i % num_width_tris != 0) { // Not a left-most triangle
				Triangle* temp = triangles[i - 1];
				t->pm1->halfedge->twin = temp->pm3->halfedge;
			}
			else {
				t->pm1->halfedge->twin = nullptr;
			}

			// Get triangle above, if it exists
			if (i >= num_width_tris) { // Not a top-most triangle
                Triangle* temp = triangles[i - num_width_tris + 1];
                t->pm3->halfedge->twin = temp->pm2->halfedge;
			}
			else {
				t->pm3->halfedge->twin = nullptr;
			}

			// Get triangle to bottom right; guaranteed to exist
			Triangle* temp = triangles[i + 1];
			t->pm2->halfedge->twin = temp->pm1->halfedge;
		}
		else {
			// Get right triangle, if it exists
			if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
				Triangle* temp = triangles[i + 1];
				t->pm3->halfedge->twin = temp->pm1->halfedge;
			}
			else {
				t->pm3->halfedge->twin = nullptr;
			}

			// Get triangle below, if it exists
			if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
				Triangle* temp = triangles[i + num_width_tris - 1];
				t->pm2->halfedge->twin = temp->pm3->halfedge;
			}
			else {
				t->pm2->halfedge->twin = nullptr;
			}

			// Get triangle to top left; guaranteed to exist
			Triangle* temp = triangles[i - 1];
			t->pm1->halfedge->twin = temp->pm2->halfedge;
		}

		topLeft = !topLeft;
	}

	clothMesh->triangles = triangles;
	this->clothMesh = clothMesh;
}
