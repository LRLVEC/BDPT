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
struct Circle
{
	vec4 plane;
	vec4 sphere;
	vec3 e1;
	uint tex;
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
layout(std430, binding = 4)buffer Circles
{
	Circle circles[];
};
void main()
{
	if (gl_GlobalInvocationID.x < circleNum)
	{
		vec3 n = normalize(circles[gl_GlobalInvocationID.x].plane.xyz);
		circles[gl_GlobalInvocationID.x].plane =
			vec4(n, -dot(n, circles[gl_GlobalInvocationID.x].sphere.xyz));
		circles[gl_GlobalInvocationID.x].e1 =
			normalize(circles[gl_GlobalInvocationID.x].e1);
	}
}