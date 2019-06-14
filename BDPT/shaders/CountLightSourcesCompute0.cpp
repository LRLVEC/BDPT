#version 450 core
#define BackwardDepth 6
#define Pi 3.14159265359
#define offset 0.005
#define originSamples 8
#define bit(a, b)			((a & (1 << uint(b))) != 0)
#define setBit(a, b)		(a |= (1 << uint(b)))
#define clearBit(a, b)		(a &= ~(1 << uint(b)))


struct LightSource
{
	uvec2 index;//geometry and num
	float a;
	float blank;
};
struct SubSource
{
	vec4 p;
	vec3 n;
	vec3 color;
};

struct Ray
{
	vec4 p0;
	vec3 n;
	float t;
};
struct Color
{
	vec3 r;
	int texR;
	vec3 t;
	int texT;
	vec3 d;
	int texD;
	vec3 g;
	int texG;
	vec3 decayFactor;
	float n;
};
struct Bound
{
	vec3 min;
	int leftChild;
	vec3 max;
	int rightChild;
};
struct BVHNode
{
	Bound bound;
	uint father;
	uint axis;
	uint geometry;
	uint geometryNum;
	vec4 blank;
};

struct TriangleGPU
{
	vec4 plane;
	vec3 p1;
	float a;
	vec3 k1;
	vec3 k2;
	vec2 uv1;
	vec2 uv2;
	vec2 uv3;
	vec2 blank;
	ivec4 nIndices;
	Color color;
};
struct Sphere
{
	vec4 sphere;
	vec3 e1;
	vec3 e2;
	Color color;
};
struct Circle
{
	vec4 plane;
	vec4 sphere;
	vec3 e1;
	int tex;
	Color color;
};
struct Cylinder
{
	vec3 c;
	float r2;
	vec3 n;
	float l;
	vec3 e1;
	Color color;
};
struct Cone
{
	vec3 c;
	float c2;
	vec3 n;
	float l2;
	vec3 e1;
	Color color;
};
struct PointLight
{
	vec3 color;
	vec3 p;
};

layout(local_size_x = 1)in;
layout(std140, row_major, binding = 1)uniform Trans
{
	mat3 trans;
	vec3 r0;
	float z0;
	float times;
};
layout(std140, binding = 3)uniform GeometryNum
{
	uint triangleNum;
	uint sphereNum;
	uint circleNum;
	uint cylinderNum;
	uint coneNum;
	uint pointLightNum;
};

layout(std430, binding = 2)buffer Triangles
{
	TriangleGPU triangles[];
};
layout(std430, binding = 3)buffer Spheres
{
	Sphere spheres[];
};
layout(std430, binding = 4)buffer Circles
{
	Circle circles[];
};
layout(std430, binding = 5)buffer Cylinders
{
	Cylinder cylinders[];
};
layout(std430, binding = 6)buffer Cones
{
	Cone cones[];
};
layout(std430, binding = 7)buffer PointLights
{
	PointLight pointLights[];
};
layout(std430, binding = 8)buffer DecayOrigin
{
	vec4 decayOrigins[originSamples];
	vec4 decayOrigin;
};
layout(std430, binding = 9)buffer BVH
{
	BVHNode bvh[];
};
layout(std430, binding = 11)buffer LightSourcesInfo
{
	uint sourcesNum;
};
layout(std430, binding = 12)buffer LightSources
{
	LightSource lightSources[];
};


void main()
{
	int c0;
	float areaAll = 0;
	uint id = 0;
	for (c0 = 0; c0 < triangleNum; ++c0)
	{
		if (triangles[c0].color.g != vec3(0))
		{
			lightSources[id].index = uvec2(2, c0);
			lightSources[id].a = areaAll +=
				triangles[c0].a;
			++id;
		}
	}
	for (c0 = 0; c0 < sphereNum; ++c0)
	{
		if (spheres[c0].color.g != vec3(0))
		{
			lightSources[id].index = uvec2(3, c0);
			lightSources[id].a = areaAll +=
				4 * Pi * spheres[c0].sphere.w;
			++id;
		}
	}
	for (c0 = 0; c0 < circleNum; ++c0)
	{
		if (circles[c0].color.g != vec3(0))
		{
			lightSources[id].index = uvec2(4, c0);
			lightSources[id].a = areaAll +=
				Pi * circles[c0].sphere.w;
			++id;
		}
	}
	for (c0 = 0; c0 < cylinderNum; ++c0)
	{
		if (cylinders[c0].color.g != vec3(0))
		{
			lightSources[id].index = uvec2(5, c0);
			lightSources[id].a = areaAll += 2 * Pi *
				sqrt(cylinders[c0].r2) * cylinders[c0].l;
			++id;
		}
	}
	for (c0 = 0; c0 < coneNum; ++c0)
	{
		if (cones[c0].color.g != vec3(0))
		{
			lightSources[id].index = uvec2(6, c0);
			lightSources[id].a = areaAll += Pi * cones[c0].l2 *
				sqrt(1 - cones[c0].c2);
			++id;
		}
	}
	sourcesNum = id;
}