
//Set the global scene parameter variables

#ifndef PARSE_VEC3_H
#define PARSE_VEC3_H
#define MAX_INPUT 256

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>


//Camera & Scene Parameters (Global Variables)
//Here we set default values, override them in parseSceneFile()

//Image Parameters
int img_width = 800, img_height = 600;
std::string imgName = "raytraced.png";

//Camera Parameters
vec3 eye = vec3(0,0,0); 
vec3 forward = vec3(0,0,-1).normalized();
vec3 up = vec3(0,1,0).normalized();
vec3 right;
float halfAngleVFOV = 35; 

//Sphere Parameters
struct Sphere {
  vec3 spherePos = vec3(0,0,2);
  float sphereRadius = 1;
  float sphere_material_index = 0;
};
Sphere spheres[MAX_INPUT];
int sphere_count = 0;

//Background Parameters
vec3 background = vec3(0,0,0);

//Material Parameters
struct Material {
  vec3 ambient = vec3(0,0,0);
  vec3 diffuse = vec3(1,1,1);
  vec3 specular = vec3(0,0,0);
  float phong_cos = 0; // the "ns"
  vec3 transmissive = vec3(0,0,0);
  float ior = 1;
};
Material materials[MAX_INPUT];
int mat_count = 0;

//Light Parameters
float point_light_r[MAX_INPUT], point_light_g[MAX_INPUT], point_light_b[MAX_INPUT];
float point_light_x[MAX_INPUT], point_light_y[MAX_INPUT], point_light_z[MAX_INPUT];
int num_lights = 0;

float dir_light_r[MAX_INPUT], dir_light_g[MAX_INPUT], dir_light_b[MAX_INPUT];
float dir_light_x[MAX_INPUT], dir_light_y[MAX_INPUT], dir_light_z[MAX_INPUT];
int num_dir_lights = 0;

vec3 ambient_light = vec3(0,0,0);
int max_depth = 5;

//store triangles

struct Triangle {
    vec3 v1, v2, v3;
    vec3 n1, n2, n3;
    bool smooth;
    int tri_material_index;
    vec3 edge1, edge2;
    vec3 N;
};

Triangle triangles[MAX_INPUT];
int tri_count = 0;

// triangle vertexes and normals
int  max_vertices;
vec3** vertices;
int max_normals;
vec3** normals;
int vert_count = 0;
int norm_count = 0;

void parseSceneFile(std::string fileName){
  //TODO: Override the default values with new data from the file "fileName"
  FILE* myFile = fopen(fileName.c_str(), "r");
  if (!myFile) {
    printf("ERROR: Could not open %s\n", fileName.c_str());
    exit(1);
  }
  char line[MAX_INPUT];
  while (fgets(line, sizeof(line), myFile)) {
    if (line[0] == '#' || line[0] == '\n') {
      continue;
    }
    if (strncmp(line, "sphere:", 7) == 0) {
      float x, y, z, r;
      if (sscanf(line + 7, "%f %f %f %f", &x, &y, &z, &r) == 4) {
        spheres[sphere_count].spherePos = vec3(x, y,z);
        spheres[sphere_count].sphereRadius = r;
        spheres[sphere_count].sphere_material_index = (mat_count > 0) ? (mat_count - 1) : 0; //avoid invalid indexing
        sphere_count++;
        printf("Sphere: (%f, %f, %f), r = %f\n", x, y, z, r);
      } else {
				std::cerr << "invalid sphere position" << std::endl;
			}
    }
    if (strncmp(line, "image_resolution:", 17) == 0) {
      int w, h;
      if (sscanf(line + 17, "%d %d", &w, &h) == 2) {
        img_width = w;
        img_height = h;
        printf("Image resloution: (%d, %d)\n", w, h);
      } else {
				std::cerr << "invalid image resolution" << std::endl;
			}
    }
    if (strncmp(line, "output_image:", 13) == 0) {
      char outputFile[MAX_INPUT];
      if (sscanf(line + 13, "%255s", outputFile) == 1) {
        imgName = std::string(outputFile);
        printf("Output image name: %s\n", imgName.c_str());
      } else {
				std::cerr << "invalid output name" << std::endl;
			}
    }
    if (strncmp(line, "camera_pos:", 11) == 0) {
      float camPos_x, camPos_y, camPos_z;
      if (sscanf(line + 11, "%f %f %f", &camPos_x, &camPos_y, &camPos_z) == 3) {
        eye = vec3(camPos_x, camPos_y, camPos_z);
        printf("Camera position: (%f, %f, %f)\n", camPos_x, camPos_y, camPos_z);
      } else {
				std::cerr << "invalid camera position" << std::endl;
			}
    }
    if (strncmp(line, "camera_fwd:", 11) == 0) {
      float fx, fy, fz;
      if (sscanf(line + 11, "%f %f %f", &fx, &fy, &fz) == 3) {
        forward = vec3(fx, fy, fz).normalized();
        printf("Camera forward direction: (%f, %f, %f)\n", fx, fy, fz);
			}
			else {
				std::cerr << "invalid camera forward direction" << std::endl;
      }
    }
    if (strncmp(line, "camera_up:", 10) == 0) {
      float ux, uy, uz;
      if (sscanf(line + 10, "%f %f %f", &ux, &uy, &uz) == 3) {
        up = vec3(ux, uy, uz).normalized();
        printf("Camera up alignment: (%f, %f, %f)\n", ux, uy, uz);
			}
			else {
				std::cerr << "invalid camera up alignment" << std::endl;
      }
    }
    if (strncmp(line, "camera_fov_ha:", 14) == 0) {
      float camFov;
      if (sscanf(line + 14, "%f", &camFov) == 1) {
        halfAngleVFOV = camFov;
        printf("Camera half-angle (vertical FOV): %f\n", camFov);
			}
			else {
				std::cerr << "invalid camera half-angle " << std::endl;
      }
    }
    if (strncmp(line, "background:", 11) == 0) {
      float backR, backG, backB;
      if (sscanf(line + 11, "%f %f %f", &backR, &backG, &backB) == 3) {
        background = vec3(backR, backG, backB);
        printf("Background color: %f %f %f\n",  backR, backG, backB);
      } else {
				std::cerr << "invalid background color" << std::endl;
      }
    }
    if (strncmp(line, "material:", 9) == 0) {
      float ar, ag, ab;
      float dr, dg, db;
      float sr, sg, sb, ns;
      float tr, tg, tb, new_ior;
      if (sscanf(line + 9, 
        "%f %f %f %f %f %f %f %f %f %f %f %f %f %f", 
        &ar, &ag, &ab, &dr, &dg, &db, &sr, &sg, &sb, &ns, &tr, &tg, &tb, &new_ior) == 14) {
      
        //storing each value into an array to deal with multiple materials in a .txt file
        materials[mat_count].ambient = vec3(ar, ag, ab);
        materials[mat_count].diffuse = vec3(dr, dg, db);
        materials[mat_count].specular = vec3(sr, sg, sb);
        materials[mat_count].transmissive = vec3(tr, tg, tb);
        materials[mat_count].phong_cos = ns;
        materials[mat_count].ior = new_ior;

        mat_count++; //up the material count
        printf("Material: %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, new_ior);
      } else {
				std::cerr << "invalid material settings" << std::endl;
      }
    }
    if (strncmp(line, "directional_light:", 18) == 0) {
      float dirR, dirG, dirB;
      float dirX, dirY, dirZ;
      if (sscanf(line + 18, "%f %f %f %f %f %f", &dirR, &dirG, &dirB, &dirX, &dirY, &dirZ) == 6) {
        // vec3 dir_light_color = vec3(dirR, dirG, dirB);
        // vec3 dir_light_direction = vec3(dirX, dirY, dirZ);
        dir_light_r[num_dir_lights] = dirR;
        dir_light_g[num_dir_lights] = dirG;
        dir_light_b[num_dir_lights] = dirB;

        dir_light_x[num_dir_lights] = dirX;
        dir_light_y[num_dir_lights] = dirY;
        dir_light_z[num_dir_lights] = dirZ;
        num_dir_lights++;
        printf("Directional light: %f %f %f %f %f %f\n", dirR, dirG, dirB, dirX, dirY, dirZ);
      } else {
				std::cerr << "invalid directional light settings" << std::endl;
      }
    }
    if (strncmp(line, "point_light:", 12) == 0) { // set position light color + pos
      float r, g, b;
      float x, y, z;
      if (sscanf(line + 12, "%f %f %f %f %f %f", &r, &g, &b, &x, &y, &z) == 6) {
        //Point light color
        point_light_r[num_lights] = r;
        point_light_g[num_lights] = g;
        point_light_b[num_lights] = b;

        //Point light position
        point_light_x[num_lights] = x;
        point_light_y[num_lights] = y;
        point_light_z[num_lights] = z;
        num_lights++;
        printf("Point light color: (%f, %f, %f)\n", r, g, b);
        printf("Point light position: (%f, %f, %f)\n", x, y, z);
      }
      else {
        std::cerr << "invalid point light settings" << std::endl;
      }
    }
    if (strncmp(line, "spot_light:", 11) == 0) {
      float r, g, b;
      float px, py, pz;
      float dx, dy, dz;
      float ang1, ang2;
      if (sscanf(line + 11, "%f %f %f %f %f %f %f %f %f %f %f",
        &r, &g, &b, &px, &py, &pz, &dx, &dy, &dz, &ang1, &ang2) == 11) {
        vec3 spotLightColor = vec3(r, g, b);
        vec3 spotLightPos = vec3(px, py, pz);
        vec3 spotLightDir = vec3(dx, dy, dz);
        float angle1 = ang1;
        float angle2 = ang2;

        printf("Spot light color: (%f, %f, %f)\n", r, g, b);
        printf("Spot light position: (%f, %f, %f)\n", px, py, pz);
        printf("Spot light direciton: (%f, %f, %f)\n", dx, dy, dz);
        printf("Angle 1: %f\n", ang1);
        printf("Angle2: %f\n", ang2);
      } else {
        std::cerr << "invalid spot light settings" << std::endl;
      }
    }
    if (strncmp(line, "ambient_light:", 14) == 0) {
      float r, g, b;
      if (sscanf(line + 14, "%f %f %f", &r, &g, &b) == 3) {
        ambient_light = vec3(r, g, b);
        printf("Ambient light color: (%f, %f, %f)\n", r, g, b);
      } else {
        std::cerr << "invalid ambient light settings" << std::endl;
      }
    }
    if (strncmp(line, "max_depth:", 10) == 0) {
      int mx_dpth;
      if (sscanf(line + 10, "%d", &mx_dpth) == 1) {
        max_depth = mx_dpth;
        printf("Max depth: %d\n", mx_dpth);
      } else {
        std::cerr << "invalid max depth" << std::endl;
      }
    }
    if (strncmp(line, "max_vertices:", 13) == 0) {
      int max_vert;
      if (sscanf(line + 13, "%d", &max_vert) == 1) {
        max_vertices = max_vert;
        vertices = new vec3*[max_vertices];
        for (int i = 0; i < max_vertices; i++) {
          vertices[i] = new vec3(0, 0, 0);
        }
        printf("Max vertices: %d\n", max_vert);
      }
      else {
        std::cerr << "invalid max vertices" << std::endl;
      }
    }
    if (strncmp(line, "max_normals:", 12) == 0) {
      int max_norm;
      if (sscanf(line + 12, "%d", &max_norm) == 1) {
        max_normals = max_norm;
        normals = new vec3*[max_normals];
        for (int i = 0; i < max_normals; i++) {
          normals[i] = new vec3(0, 0, 0);
        }
        printf("Max normals: %d\n", max_norm);
      }
      else {
        std::cerr << "invalid max normals" << std::endl;
      }
    }
    if (strncmp(line, "vertex:", 7) == 0) {
      if (max_vertices == 0) {
        std::cerr << "must specify max vertices" << std::endl;
      }
      int x, y, z;
      if (sscanf(line + 7, "%d %d %d", &x, &y, &z) == 3) {
        vec3* vertex = new vec3(x, y, z);
        vertices[vert_count] = vertex;
        vert_count++;
        printf("Vertex: (%d, %d, %d)\n", x, y, z);
      }
      else {
        std::cerr << "invalid vertex" << std::endl;
      }
    }
    if (strncmp(line, "normal:", 7) == 0) {
      if (max_normals == 0) {
        std::cerr << "must specify max normals" << std::endl;
      }
      int x, y, z;
      if (sscanf(line + 7, "%d %d %d", &x, &y, &z) == 3) {
        vec3* normal = new vec3(x, y, z);
        normals[norm_count] = normal;
        norm_count++;
        printf("Normal: (%d, %d, %d)\n", x, y, z);
      }
      else {
        std::cerr << "invalid normal" << std::endl;
      }
    }
    if (strncmp(line, "triangle:", 9) == 0) {
      int v1, v2, v3;
      if (sscanf(line + 9, "%d %d %d", &v1, &v2, &v3) == 3) {
        triangles[tri_count].v1 = *vertices[v1];
        triangles[tri_count].v2 = *vertices[v2];
        triangles[tri_count].v3 = *vertices[v3];

        triangles[tri_count].edge1 = triangles[tri_count].v2 - triangles[tri_count].v1;
        triangles[tri_count].edge2 = triangles[tri_count].v3 - triangles[tri_count].v1;
        triangles[tri_count].N = cross(triangles[tri_count].edge1, triangles[tri_count].edge2).normalized();
        triangles[tri_count].smooth = false;
        triangles[tri_count].tri_material_index = (mat_count > 0) ? (mat_count - 1) : 0;
        tri_count++;
        printf("Triangle vertex numbers: %d %d %d\n", v1, v2, v3);
      }
      else {
        std::cerr << "invalid triangle" << std::endl;
      }
    }
    if (strncmp(line, "normal_triangle:", 16) == 0){
      int v1, v2, v3;
      int n1, n2, n3;
      if (sscanf(line + 16, "%d %d %d %d %d %d", &v1, &v2, &v3, &n1, &n2, &n3) == 6) {
        triangles[tri_count].v1 = *vertices[v1];
        triangles[tri_count].v2 = *vertices[v2];
        triangles[tri_count].v3 = *vertices[v3];
        triangles[tri_count].n1 = *normals[n1];
        triangles[tri_count].n2 = *normals[n2];
        triangles[tri_count].n3 = *normals[n3];
        triangles[tri_count].smooth = true;
        tri_count++;
        printf("Normal triangle vertex and normal numbers: %d %d %d %d %d %d\n", v1, v2, v3, n1, n2, n3);
      }
      else {
        std::cerr << "invalid normal triangle" << std::endl;
      }
    }
  }
  printf("Orthogonal Camera Basis:\n");
  forward = forward.normalized();
  right = cross(up, forward).normalized();
  up = cross(forward, right).normalized();

  printf("forward: %f,%f,%f\n",forward.x,forward.y,forward.z);
  printf("right: %f,%f,%f\n",right.x,right.y,right.z);
  printf("up: %f,%f,%f\n",up.x,up.y,up.z);

  fclose(myFile);
}

#endif