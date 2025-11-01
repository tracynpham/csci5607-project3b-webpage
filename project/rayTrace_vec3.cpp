
//To Compile: g++ -fsanitize=address -std=c++14 rayTrace_vec3.cpp -o ray
//Compile with Parallelism on MacOS: g++-15 -fopenmp -std=c++14 rayTrace_vec3.cpp -o ray
// v2 for csci 5607 project 3b
// tracy and mary 

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS // For fopen and sscanf
#define _USE_MATH_DEFINES 
#endif

//Images Lib includes:
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "image_lib.h" //Defines an image class and a color class

//#Vec3 Library
#include "vec3.h"

//High resolution timer
#include <chrono>
#include <omp.h>

//Scene file parser
#include "parse_vec3.h"
#include <algorithm>

//added struct to hold information about the ray intersection with the sphere
struct HitInformation {
  vec3 point; //hit point
  vec3 normal; //normal along hit point
  int sphere_num = -1; //index of the hit sphere
  int tri_num = -1;
  bool hit = false;
  float dist = -1.0f;
};


Color ApplyLightingModel(vec3 start, vec3 dir, HitInformation& hitInfo, int depth = 0);
Color evaluateRayTree(vec3 start, vec3 dir, int depth = 0);

//returns were exactly the ray intersected the sphere
float whereRaySphereIntersect(vec3 start, vec3 dir, vec3 center, float radius){
  float a = dot(dir,dir); 
  vec3 toStart = (start - center);
  float b = 2 * dot(dir,toStart);
  float c = dot(toStart,toStart) - radius*radius;
  float discr = b*b - 4*a*c;
  if (discr < 0) return -1;
  else{
    float t0 = (-b + sqrt(discr))/(2*a);
    float t1 = (-b - sqrt(discr))/(2*a);
    if (t0 > 0 && t1 > 0) 
      return (t0 < t1) ? t0 : t1; //only return the smaller distance
    else if (t0 > 0) //t1 was negative 
      return t0;
    else if (t1 > 0) //t0 was negative
      return t1;
    else
      return -1; //no intersection
  }
}
bool sameSide(vec3 p1, vec3 p2, vec3 a, vec3 b) {
  vec3 cp1 = cross(b-a, p1-a);
  vec3 cp2 = cross(b-a, p2-a);
  return dot(cp1, cp2) >= 0;
}

bool pointInRegTriangle(vec3 p, vec3 a, vec3 b, vec3 c) {
  return sameSide(p, a, b, c) && sameSide(p, b, a, c) && sameSide(p, c, a, b);
}

// using cramer's rule from the textbook Real-Time Collision Detection
// https://ceng2.ktu.edu.tr/~cakir/files/grafikler/Texture_Mapping.pdf
bool pointInNormalTriangle(vec3 p, vec3 a, vec3 b, vec3 c, float &alpha, float &beta, float &gamma) {
  vec3 v0 = b - a;
  vec3 v1 = c - a;
  vec3 v2 = p - a;

  float d00 = dot(v0, v0);
  float d01 = dot(v0, v1);
  float d11 = dot(v1, v1);
  float d20 = dot(v2, v0);
  float d21 = dot(v2, v1);
  float denom = (d00 * d11) - (d01 * d01);

  float EPSILON = 1e-6f;
  if (fabs(denom) < EPSILON) {
    alpha = beta = gamma = -1;
    return false;
  }

  beta = (d11 * d20 - d01 * d21) / denom;
  gamma = (d00 * d21 - d01 * d20) / denom;
  alpha = 1.0f - beta - gamma;

  return (alpha >= 0.0f) && (beta >= 0.0f) && (gamma >= 0.0f);
}

void TriangleIntersection(vec3 start, vec3 dir, HitInformation& hitInfo, float& closest_dist) {
  for (int tri = 0; tri < tri_count; tri++) {
    vec3 triangle_v1 = triangles[tri].v1;
    vec3 triangle_v2 = triangles[tri].v2;
    vec3 triangle_v3 = triangles[tri].v3;
    // calculate point first
    vec3 c0 = triangle_v1;
    vec3 N;
    vec3 point;
    // depending on type triangle, after we find the point, then calculate its normal
    // this is for regular triangles
    if (!triangles[tri].smooth) {
      vec3 edge1 = triangles[tri].edge1;
      vec3 edge2 = triangles[tri].edge2;
      N = triangles[tri].N;
      //ray plane intersection formula
      if (fabs(dot(dir, N)) < 0.001f) { //checks if ray is parallel to plane // TEMP CHANGE
        continue;
      }
      float d = -dot(c0, N);
      float t = -(dot(start, N) + d) / dot(dir, N);
      if (t < 0) { //intersection is behind the camera
        continue; //ignore
      } 
        
      point = start+t*dir; //point where ray intersect plane
      if (pointInRegTriangle(point, triangle_v1, triangle_v2, triangle_v3)) {
        if (closest_dist < 0 || t < closest_dist) {
          closest_dist = t;
          hitInfo.tri_num = tri;
          hitInfo.hit = true;
          hitInfo.point = point;
          if (dot(dir, N) > 0) { N = N * -1; }
          hitInfo.normal = N;
          hitInfo.dist = t;
        }
      } 
    }
    // this is for normal triangles
    else {
      vec3 edge1 = triangles[tri].edge1;
      vec3 edge2 = triangles[tri].edge2;
      //N = triangles[tri].N;
      // TEMP CHANGE
      vec3 planeNormal = triangles[tri].N.normalized();
      //ray plane intersection formula
      //if (fabs(dot(dir, N)) < 0.001f) { //checks if ray is parallel to plane 
      //  continue;
      //}
      float d = -dot(c0, planeNormal);
      float t = -(dot(start, planeNormal) + d) / dot(dir, planeNormal);
      if (t < 0) { //intersection is behind the camera
        continue; //ignore
      }
      point = start + t * dir; //point where ray intersect plane
      float alpha, beta, gamma;
      if (pointInNormalTriangle(point, triangle_v1, triangle_v2, triangle_v3, alpha, beta, gamma)) {
        if (closest_dist < 0 || t < closest_dist) {
          closest_dist = t;
          hitInfo.tri_num = tri;
          hitInfo.hit = true;
          hitInfo.point = point;

          vec3 n1 = triangles[tri].n1.normalized();
          vec3 n2 = triangles[tri].n2.normalized();
          vec3 n3 = triangles[tri].n3.normalized();
          vec3 interpolatedNormal = (alpha * n1 + beta * n2 + gamma * n3).normalized();
          // if (dot(dir, interpolatedNormal) > 0) {
          //   interpolatedNormal = interpolatedNormal * -1;
          // }
          //TEMP CHANGE
          if (dot(planeNormal, interpolatedNormal) < 0){
            interpolatedNormal = interpolatedNormal * -1;
          }
          hitInfo.normal = interpolatedNormal;
          hitInfo.dist = t;
        }
      }
    }

  }
}

//Find intersection of different shapes and save into hitInfo
bool FindIntersection(vec3 start, vec3 dir, HitInformation& hitInfo) {
  float closest_dist = -1;
  for (int s = 0; s < sphere_count; s++) {
    vec3 spherePos = spheres[s].spherePos;
    float radius = spheres[s].sphereRadius;
    float dist = whereRaySphereIntersect(start, dir, spherePos, radius);
    if (dist > 0 && (closest_dist < 0 || dist < closest_dist)) {
      closest_dist = dist;
      hitInfo.sphere_num = s;
      hitInfo.hit = true;
      hitInfo.point = start + dir * dist;
      hitInfo.normal = (hitInfo.point - spherePos).normalized();
      hitInfo.dist = dist;
    }
  }
  TriangleIntersection(start, dir, hitInfo, closest_dist);
  return hitInfo.hit;
}
Color evaluateRayTree(vec3 start, vec3 dir, int depth) {
  //base case
  if (depth > max_depth) {
    return Color(0,0,0);
  }
  bool hit_something = false;
  HitInformation hit;
  hit_something = FindIntersection(start, dir, hit);
  if (hit_something) {
    return ApplyLightingModel(start, dir, hit, depth);
  } else {
    return Color(background.x, background.y, background.z);
  }
}

vec3 reflect(vec3 d, vec3 n) {
  return (d-2.0f*dot(d, n)*n);
}

bool refract(vec3 d, vec3 n, float r, vec3& t) {
  d = d.normalized();
  n = n.normalized();
  float cos_theta = -dot(n, d);
  float k = 1.0f - r * r * (1.0f - cos_theta * cos_theta);
  if (k < 0.0f) return false; 
  t = (r * d) + ((r * cos_theta) - std::sqrt(k)) * n;
  return true;
}


Color ApplyLightingModel(vec3 start, vec3 dir, HitInformation& hitInfo, int depth) {
  vec3 contribution = vec3(0, 0, 0);
  int mat_index;
  if (hitInfo.tri_num >= 0)
    mat_index = triangles[hitInfo.tri_num].tri_material_index;
  else
    mat_index = spheres[hitInfo.sphere_num].sphere_material_index;
  // getting diffuse, specular, ambient, and transmissive coefficients from parser
  vec3 kd = materials[mat_index].diffuse;
  vec3 ks = materials[mat_index].specular;
  vec3 ka = materials[mat_index].ambient;
  vec3 kt = materials[mat_index].transmissive;
  float ns = materials[mat_index].phong_cos;

  // adding the contribution from the ambient lighting once , no need to loop through the lights to get it 
  contribution = contribution + (ka * ambient_light);

  // compute view vector pointing from hitpoint to the camera
  vec3 V = (eye - hitInfo.point).normalized();

  //POINT LIGHT
  for (int i = 0; i < num_lights; i++) {
    // get light position and color, compute light vector pointing from hitpoint towards the light 
    vec3 lightPos = vec3(point_light_x[i], point_light_y[i], point_light_z[i]);
    vec3 lightColor = vec3(point_light_r[i], point_light_g[i], point_light_b[i]);
    vec3 L = (lightPos - hitInfo.point).normalized();

    // shadow ray, small offset to avoid shadows intersectiong the surface
    vec3 shadow = hitInfo.point + hitInfo.normal * 0.001f;
    HitInformation shadow_hit;
    // check if light is blocked by an object
    bool blocked = FindIntersection(shadow, L, shadow_hit);
    float light_distance = (lightPos - hitInfo.point).length(); 
    if (blocked && (shadow_hit.point - shadow).length() < light_distance){
      // this means light was blocked, skip adding diffuse and specular contributions from this light
      continue;
    }
    
    // diffuse lighting
    float diff = std::max(0.0f, dot(hitInfo.normal, L)); // n dot L
    vec3 diffuse = kd * diff;

    //specular lighting (blinn-phong)
    vec3 H = (L + V).normalized();
    float spec_angle = std::max(0.0f, dot(hitInfo.normal, H));
    vec3 specular = ks * pow(spec_angle, ns);

    // combine diffuse + specular contributions
    float attenuation = 1.0f / (light_distance * light_distance);
    vec3 totalLight = (diffuse + specular) * lightColor * attenuation;

    // add this light's contribution to the running total
    contribution = contribution + totalLight * (1.0f / num_lights);
  }
  //DIRECTIONAL LIGHT
  for (int i = 0; i < num_dir_lights; i++) {
    vec3 L = vec3(dir_light_x[i], dir_light_y[i], dir_light_z[i]).normalized() * -1;
    vec3 lightColor = vec3(dir_light_r[i], dir_light_g[i], dir_light_b[i]);

    vec3 shadow = hitInfo.point + hitInfo.normal * 0.001f;
    HitInformation shadow_hit;
    // check if light is blocked by an object
    bool blocked = FindIntersection(shadow, L, shadow_hit);
    if (blocked) {
      continue;
    }
    // diffuse lighting
    float diff = std::max(0.0f, dot(hitInfo.normal, L)); // n dot L
    vec3 diffuse = kd * diff;

    //specular lighting (blinn-phong)
    vec3 H = (L + V).normalized();
    float spec_angle = std::max(0.0f, dot(hitInfo.normal, H));
    vec3 specular = ks * pow(spec_angle, ns);

    vec3 totalLight = (diffuse + specular) * lightColor;

    // add this light's contribution to the running total
    contribution = contribution + totalLight;
  }
  //SPOT LIGHT
  for (int i = 0; i < spot_count; i++) {
    vec3 spotColor = spotLights[i].color; 
    vec3 spotPos = spotLights[i].pos;      
    vec3 spotDir = spotLights[i].dir;      
    float ang1 = cosf(spotLights[i].angle1);
    float ang2 = cosf(spotLights[i].angle2);

    vec3 L = spotPos - hitInfo.point;
    float distToLightSqr = dot(L, L);
    vec3 LNorm = L.normalized();
    float spotAngle = dot(LNorm * -1, spotDir); 

    float spotIntensity = 1.0f / distToLightSqr;
    if (spotAngle > ang1) {
      spotIntensity *= 1.0f; //acts light point light
    }
    else if (spotAngle > ang2) { //between both angles for smooth fall off
      float falloffRange = ang1 - ang2;
      spotIntensity *= (spotAngle - ang2) / falloffRange;
    }
    else { 
      spotIntensity = 0.0f; //no contribution
    }
    //dealing with shadows
    HitInformation shadowHit;
    vec3 shadowOrigin = hitInfo.point + hitInfo.normal * 0.001f;
    bool blocked = FindIntersection(shadowOrigin, LNorm, shadowHit);
    float lightDistance = sqrt(distToLightSqr);
    if (!blocked || shadowHit.dist >= lightDistance) {
      //diffuse
      float diff = std::max(0.0f, dot(hitInfo.normal, LNorm));
      vec3 diffuse = kd * diff;

      //Specular (Blinn-Phong)
      vec3 V = (eye - hitInfo.point).normalized();
      vec3 H = (LNorm + V).normalized();
      float specAngle = std::max(0.0f, dot(hitInfo.normal, H));
      vec3 specular = ks * pow(specAngle, ns);

      // Add contribution
      contribution = contribution + (diffuse + specular) * spotColor * spotIntensity;
    }
  }
  //Following the pseudocode from lecture slides 13
  if (depth <= max_depth) {
    float mat_ior = materials[mat_index].ior;
    bool entering = dot(dir, hitInfo.normal) < 0.0f; // determine if ray is entering or exiting the object
      
    // set surface normal to point in correct direction for refractio
    // flip it when exiting 
    vec3 n = hitInfo.normal;
    vec3 nRefract = entering ? n : vec3(-n.x, -n.y, -n.z);

    // set refractive indices for snell's law
    float n1 = entering ? 1.0f : mat_ior;
    float n2 = entering ? mat_ior : 1.0f;
    float snell_ratio = n1 / n2;

    // reflect ray computation (Fresnell effect)
      vec3 r = reflect(dir, nRefract).normalized();
      vec3 reflect_start = hitInfo.point + nRefract * 0.001f;
      Color reflect_color = evaluateRayTree(reflect_start, r, depth + 1);
      vec3 reflect_contrib = vec3(reflect_color.r, reflect_color.g, reflect_color.b) * ks;

      // refraction ray (only if refract is true)
      vec3 t;
      bool canRefract = refract(dir, nRefract, snell_ratio, t); //snell's law
      vec3 k = vec3(0.0f, 0.0f, 0.0f);
      if (canRefract) {
        vec3 t_dir = t.normalized();
        vec3 refract_start = hitInfo.point + t_dir * 0.001f;
        Color refract_color = evaluateRayTree(refract_start, t.normalized(), depth + 1);
        k = vec3(refract_color.r, refract_color.g, refract_color.b) * kt;
      }
      
      contribution = contribution + (reflect_contrib + (1.0f) * (k));
    }
  contribution.clampTo1(); // clamp so none of the exponents exceed 1
  return Color(contribution.x, contribution.y, contribution.z); // this is where i converted it to a color
}

int main(int argc, char** argv){

  //Read command line paramaters to get scene file
  if (argc != 2){
    std::cout << "Usage: ./a.out scenefile\n";
    return(0);
  }
  std::string secenFileName = argv[1];

  //Parse Scene File
  parseSceneFile(secenFileName);

  float imgW = img_width, imgH = img_height;
  float halfW = imgW/2, halfH = imgH/2;
  float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));
  Image outputImgSerial = Image(img_width,img_height);
  
  std::cout << "--------------------Performance-------------------- " << std::endl;
 /* auto t_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < img_width; i++) {
      for (int j = 0; j < img_height; j++) {
        float u = (halfW - (imgW) * ((i + 0.5f) / imgW));
        float v = (halfH - (imgH) * ((j + 0.5f) / imgH));
        vec3 p = eye - d * forward + u * right + v * up;
        vec3 rayDir = (p - eye).normalized();
        Color color = evaluateRayTree(eye, rayDir, 0);
        outputImgSerial.setPixel(i, j, color);
      }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  double serial_duration = std::chrono::duration<double>(t_end - t_start).count();
  printf("Serial Rendering took %.2f s\n", serial_duration);
  outputImgSerial.write(imgName.c_str());*/

  //parallel render using OpenMP
  Image outputImgParallel = Image(img_width,img_height);
  auto t_start_parallel = std::chrono::high_resolution_clock::now();
  #pragma omp parallel for
  for (int i = 0; i < img_width; i++) {
    for (int j = 0; j < img_height; j++) {
      float u = (halfW - (imgW) * ((i + 0.5f) / imgW));
      float v = (halfH - (imgH) * ((j + 0.5f) / imgH));
      vec3 p = eye - d * forward + u * right + v * up;
      vec3 rayDir = (p - eye).normalized();
      Color color = evaluateRayTree(eye, rayDir, 0);
      outputImgParallel.setPixel(i, j, color);
      }
  }
  auto t_end_parallel = std::chrono::high_resolution_clock::now();
  double parallel_duration = std::chrono::duration<double>(t_end_parallel - t_start_parallel).count();
  printf("Parallel Rendering took %.2f s\n", parallel_duration);
  outputImgParallel.write(imgName.c_str());

  //double speedup = serial_duration / parallel_duration;
  //std::cout << "Speedup: " << speedup << std::endl;
  std::cout << "----------End of Performance Analysis---------- " << std::endl;

  // check if vertices and/or normals exist
  // if they do, delete them to free up the space we allocate
  if (vertices) {
    for (int i = 0; i < max_vertices; i++) {
      delete vertices[i];
    }
    delete[] vertices;
    printf("Freed up space from vertices array\n");
  }
  if (normals) {
    for (int i = 0; i < max_normals; i++) {
      delete normals[i];
    }
    delete[] normals;
    printf("Freed up space from normals array\n");
  }
  return 0;
}