#version 450 core
layout(local_size_x = 1024)in;

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
	vec3 blank;
	float n;
};
struct Triangle
{
	vec3 p1;
	vec3 p2;
	vec3 p3;
	vec2 uv1;
	vec2 uv2;
	vec2 uv3;
	vec2 blank;
	ivec4 nIndices;
	Color color;
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
layout(std140, binding = 3)uniform GeometryNum
{
	uint triangleNum;
	uint sphereNum;
	uint circleNum;
	uint cylinderNum;
	uint coneNum;
	uint pointLightNum;
};

layout(std430, binding = 1)buffer TriangleOrigin
{
	Triangle trianglesOrigin[];
};
layout(std430, binding = 2)buffer Triangles
{
	TriangleGPU triangles[];
};

void main()
{
	if (gl_GlobalInvocationID.x < triangleNum)
	{
		triangles[gl_GlobalInvocationID.x].p1 = trianglesOrigin[gl_GlobalInvocationID.x].p1;
		vec3 e1 =
			trianglesOrigin[gl_GlobalInvocationID.x].p2 -
			trianglesOrigin[gl_GlobalInvocationID.x].p1;
		vec3 e2 =
			trianglesOrigin[gl_GlobalInvocationID.x].p3 -
			trianglesOrigin[gl_GlobalInvocationID.x].p1;
		vec3 n = cross(e1, e2);
		float s = dot(n, n);
		n = normalize(n);
		triangles[gl_GlobalInvocationID.x].plane =
			vec4(n, -dot(n, triangles[gl_GlobalInvocationID.x].p1));
		float d = dot(e1, e2);
		triangles[gl_GlobalInvocationID.x].a = sqrt(s) / 2;
		triangles[gl_GlobalInvocationID.x].k1 = (dot(e2, e2) * e1 - d * e2) / s;
		triangles[gl_GlobalInvocationID.x].k2 = (dot(e1, e1) * e2 - d * e1) / s;
		triangles[gl_GlobalInvocationID.x].uv1 = trianglesOrigin[gl_GlobalInvocationID.x].uv1;
		triangles[gl_GlobalInvocationID.x].uv2 = trianglesOrigin[gl_GlobalInvocationID.x].uv2;
		triangles[gl_GlobalInvocationID.x].uv3 = trianglesOrigin[gl_GlobalInvocationID.x].uv3;
		triangles[gl_GlobalInvocationID.x].nIndices = trianglesOrigin[gl_GlobalInvocationID.x].nIndices;
		triangles[gl_GlobalInvocationID.x].color = trianglesOrigin[gl_GlobalInvocationID.x].color;
	}
}