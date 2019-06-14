#version 450 core
layout(local_size_x = 1)in;
#define offset 0.001
#define originSamples 8
#define bit(a, b)			((a & (1 << uint(b))) != 0)
#define setBit(a, b)		(a |= (1 << uint(b)))
#define clearBit(a, b)		(a &= ~(1 << uint(b)))

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
layout(std430, binding = 8)buffer DecayOrigin
{
	vec4 decayOrigins[originSamples];
	vec4 decayOrigin;
};
layout(std430, binding = 9)buffer BVH
{
	BVHNode bvh[];
};

float getPlaneT(Ray ray, vec4 plane)
{
	return -dot(plane, ray.p0) / dot(plane.xyz, ray.n);
}
vec2 getTriangleUV(vec3 pos, uint num)
{

	vec3 d = pos - triangles[num].p1;
	return vec2
	(
		dot(d, triangles[num].k1),
		dot(d, triangles[num].k2)
	);
}
bool triangleTest(vec2 uv)
{
	if (all(greaterThanEqual(uv, vec2(0, 0))) && (uv.x + uv.y <= 1.0f))
		return true;
	return false;
}
bool judgeHitBox(Ray ray, Bound bound)
{
	if (all(lessThanEqual(bound.min, ray.p0.xyz)) && all(lessThanEqual(ray.p0.xyz, bound.max)))
		return true;
	vec3 tmin = (bound.min - ray.p0.xyz) / ray.n;
	vec3 tmax = (bound.max - ray.p0.xyz) / ray.n;
	vec3 mintt = min(tmin, tmax);
	tmax = max(tmin, tmax);
	float maxt = max(mintt.x, max(mintt.y, mintt.z));
	if (maxt > 0)
	{
		float mint = min(tmax.x, min(tmax.y, tmax.z));
		return (maxt <= mint) && (ray.t < 0 || maxt <= ray.t);
	}
	else
		return false;
}
float random(vec2 st)
{
	return 2 * fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123) - 1;
}


void main()
{
	Ray ray;
	ray.p0 = vec4(r0, 1);
	ray.n = normalize
	(vec3
	(
		random(r0.xy + vec2(gl_GlobalInvocationID.x & 1)),
		random(r0.yz + vec2((gl_GlobalInvocationID.x & 2) / 2)),
		random(r0.zx + vec2((gl_GlobalInvocationID.x & 4) / 4))
	));
	vec3 decayTemp = vec3(0);
	float nTemp = 1;
	while (true)
	{
		float t = -1;
		uint n = 0;
		vec3 decayNow = vec3(0);
		float nNow;
		ray.t = -1;
		uint now = 0;
		uint bvhStack = 0;
		uint bvhState;
		int bvhSP = -1;
		for (;;)
		{
			if (judgeHitBox(ray, bvh[now].bound))
			{
				uint geometry = bvh[now].geometry;
				if (geometry != 0)
				{
					uint n = bvh[now].geometryNum;
					switch (geometry)
					{
						/*	for (n = 0; n < planeNum; ++n)
								{
									float tt = getPlaneT(ray, planes[n].plane);
									if (tt > 0 && (tt < ray.t || ray.t < 0))
									{
										ray.t = tt;
										vec3 p1 = ray.p0.xyz + ray.n * ray.t;
										tempColor = planes[n].color;
										tempColor.g = ((int(4.2 * p1.x) + int(4.2 * p1.y)) % 2u) * tempColor.g;
										tempN = planes[n].plane.xyz;
									}
							}*/
						case 2:
						{
							float tt = getPlaneT(ray, triangles[n].plane);
							if (tt > 0 && (tt < ray.t || ray.t < 0))
							{
								vec2 uv = getTriangleUV(ray.p0.xyz + ray.n * tt, n);
								if (triangleTest(uv))
								{
									ray.t = tt;
									float dn = dot(ray.n, triangles[n].plane.xyz);
									decayNow = triangles[n].color.decayFactor * sign(dn);
									nNow = dn > 0 ? triangles[n].color.n : 1 / triangles[n].color.n;
								}
							}
							break;
						}
						case 3:
						{
							vec3 d = spheres[n].sphere.xyz - ray.p0.xyz;
							float s = spheres[n].sphere.w - dot(cross(d, ray.n), cross(d, ray.n));
							if (s >= 0)
							{
								s = sqrt(s);
								float k = dot(d, ray.n);
								float tt = -1;
								if (k + s > 0)tt = k + s;
								if (k - s > 0)tt = k - s;
								if (tt > 0 && (tt < ray.t || ray.t < 0))
								{
									ray.t = tt;
									float dn = dot(ray.n, (ray.p0.xyz + ray.t * ray.n - spheres[n].sphere.xyz) / sqrt(spheres[n].sphere.w));
									decayNow = spheres[n].color.decayFactor * sign(dn);
									nNow = dn > 0 ? spheres[n].color.n : 1 / spheres[n].color.n;
								}
							}
							break;
						}
						case 4:
						{
							float tt = getPlaneT(ray, circles[n].plane);
							if (tt > 0 && (tt < ray.t || ray.t < 0))
							{
								vec3 d = ray.p0.xyz + ray.n * tt - circles[n].sphere.xyz;
								if (dot(d, d) <= circles[n].sphere.w)
								{
									ray.t = tt;
									float dn = dot(ray.n, circles[n].plane.xyz);
									decayNow = circles[n].color.decayFactor * sign(dn);
									nNow = dn > 0 ? circles[n].color.n : 1 / circles[n].color.n;
								}
							}
							break;
						}
						case 5:
						{
							float nn0 = dot(ray.n, cylinders[n].n);
							float cnn02 = 1 - nn0 * nn0;
							if (cnn02 != 0)
							{
								vec3 d = ray.p0.xyz - cylinders[n].c;
								float nd = dot(cylinders[n].n, d);
								vec3 j = d - nd * cylinders[n].n;
								float n0j = dot(ray.n, j);
								float k = n0j * n0j + cnn02 * (cylinders[n].r2 - dot(j, j));
								if (k > 0)
								{
									k = sqrt(k);
									float tt = -1;
									float v;
									if (k - n0j > 0)
									{
										tt = (k - n0j) / cnn02;
										v = nd + nn0 * tt;
										if (v > cylinders[n].l || v < 0)
											tt = -1;
									}
									if (k + n0j < 0)
									{
										float ttt = -(k + n0j) / cnn02;
										float ut = nd + nn0 * ttt;
										if (ut <= cylinders[n].l && ut >= 0)
										{
											tt = ttt;
											v = ut;
										}
									}
									if (tt > 0 && (tt < ray.t || ray.t < 0))
									{
										ray.t = tt;
										float dn = dot(ray.n, normalize(d + ray.n * ray.t - cylinders[n].n * v));
										decayNow = cylinders[n].color.decayFactor * sign(dn);
										nNow = dn > 0 ? cylinders[n].color.n : 1 / cylinders[n].color.n;
									}
								}
							}
							break;
						}
						case 6:
						{
							vec3 d = ray.p0.xyz - cones[n].c;
							float nn0 = dot(ray.n, cones[n].n);
							float dn0 = dot(d, cones[n].n);
							float dn = dot(d, ray.n);
							float d2 = dot(d, d);
							float a = cones[n].c2 - nn0 * nn0;
							float b = nn0 * dn0 - dn * cones[n].c2;
							float c = d2 * cones[n].c2 - dn0 * dn0;
							float s = b * b - a * c;
							if (s > 0)
							{
								s = sqrt(s);
								float tt = -1;
								float r2;
								tt = (b + s) / a;
								if (tt > 0)
								{
									r2 = d2 + tt * tt + 2 * dn * tt;
									float k = dn0 + nn0 * tt;
									if (r2 > cones[n].l2 || k < 0)tt = -1;
								}
								float ttt = (b - s) / a;
								if (ttt > 0 && (ttt < tt || (tt < 0)))
								{
									float r2t = d2 + ttt * ttt + 2 * dn * ttt;
									float k = dn0 + nn0 * ttt;
									if (r2t <= cones[n].l2 && k > 0)
									{
										tt = ttt;
										r2 = r2t;
									}
								}
								if (tt > 0 && (tt < ray.t || ray.t < 0))
								{
									ray.t = tt;
									float dn = dot(ray.n, normalize(d + ray.n * ray.t - cones[n].n * sqrt(r2 / cones[n].c2)));
									decayNow = cones[n].color.decayFactor * sign(dn);
									nNow = dn > 0 ? cones[n].color.n : 1 / cones[n].color.n;
								}
							}
							break;
						}
					}
				}
				if (bvh[now].bound.leftChild != 0)
				{
					setBit(bvhStack, ++bvhSP);
					if (ray.n[bvh[now].axis] >= 0 || bvh[now].bound.rightChild == 0)
					{
						++now;
						clearBit(bvhState, bvhSP);
					}
					else
					{
						now = bvh[now].bound.rightChild;
						setBit(bvhState, bvhSP);
					}
					continue;
				}
			}
			while (bvhSP >= 0)
			{
				now = bvh[now].father;
				if (bit(bvhStack, bvhSP))
				{
					if (bit(bvhState, bvhSP))
					{
						clearBit(bvhStack, bvhSP);
						++now;
						break;
					}
					else if (bvh[now].bound.rightChild != 0)
					{
						clearBit(bvhStack, bvhSP);
						now = bvh[now].bound.rightChild;
						break;
					}
				}
				--bvhSP;
			}
			if (bvhSP < 0)break;
		}
		if (ray.t < 0)
		{
			decayOrigins[gl_GlobalInvocationID.x] = vec4(decayTemp, nTemp);
			return;
		}
		else
		{
			decayTemp += decayNow;
			nTemp *= nNow;
			ray.p0 += vec4(ray.n * (ray.t + offset), 0);
		}
	}
}