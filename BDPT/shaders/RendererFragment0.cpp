#version 450 core
#define BackwardDepth 5
#define ForwardDepth 5
#define SubSoureMax 5
#define Pi 3.14159265359
#define offset 0.0001
#define minColor 0.01
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
struct HitInfo
{
	vec3 n;
	float t;
	vec2 uv;
	uvec2 hitObj;
};
struct ForwardStack
{
	vec4 tp0;
	vec3 tn;
	int t;
	vec3 tColor;
	vec3 tDecayFactor;
	vec4 dp0;
	vec3 dn;
	int d;
	vec3 dColor;
	vec3 dDecayFactor;
	int depth;
};
struct BackwardStack
{
	vec4 tp0;
	vec3 tn;
	int t;
	vec3 tRatio;
	vec3 tDecayFactor;
	vec4 dp0;
	vec3 dn;
	int d;
	vec3 dRatio;
	vec3 dDecayFactor;
	int depth;
};

//to do:	check about the decay in light sources...
//			why it is not correct...

layout(std140, binding = 0)uniform Size
{
	uvec2 size;
};
layout(std140, row_major, binding = 1)uniform Trans
{
	mat3 trans;
	vec3 r0;
	float z0;
	float times;
};
layout(binding = 2, rgba32f)uniform image2D image;
layout(binding = 1)uniform sampler2DArray texSmp;
layout(binding = 2)uniform samplerCube cubeSmp;
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
layout(std430, binding = 10)buffer Counter
{
	uint counter;
};
layout(std430, binding = 11)buffer LightSourcesInfo
{
	uint sourcesNum;
};
layout(std430, binding = 12)buffer LightSources
{
	LightSource lightSources[];
};

SubSource subSources[SubSoureMax];
int subID;
float sa;

float random(vec2 st)
{
	return  fract(sin(dot(st.xy, vec2(12.9898, 78.233))) * 43758.5453123);
}
vec3 getLambertDirection(vec3 n, vec2 rdUV)
{
	vec3 u;
	vec3 v;
	rdUV.y *= 2.0 * Pi;
	if (abs(n.x) > 0.7)
	{
		float s = sqrt(1 - n.y * n.y);
		u = vec3(-n.z, 0, n.x) / s;
		v = vec3(n.x * n.y / s, -s, n.y * n.z / s);
	}
	else
	{
		float s = sqrt(1 - n.x * n.x);
		u = vec3(0, n.z, -n.y) / s;
		v = vec3(-s, n.x * n.y / s, n.x * n.z / s);
	}
	float sinTheta = sqrt(1 - rdUV.x);
	return sinTheta * cos(rdUV.y) * u + sinTheta * sin(rdUV.y) * v + sqrt(rdUV.x) * n;
}
vec3 getDirection(vec3 n, vec2 rdUV)
{
	vec3 u;
	vec3 v;
	rdUV.y *= 2.0 * Pi;
	if (abs(n.x) > 0.7)
	{
		float s = sqrt(1 - n.y * n.y);
		u = vec3(-n.z, 0, n.x) / s;
		v = vec3(n.x * n.y / s, -s, n.y * n.z / s);
	}
	else
	{
		float s = sqrt(1 - n.x * n.x);
		u = vec3(0, n.z, -n.y) / s;
		v = vec3(-s, n.x * n.y / s, n.x * n.z / s);
	}
	float sinTheta = sqrt(1 - rdUV.x * rdUV.x);
	return sinTheta * cos(rdUV.y) * u + sinTheta * sin(rdUV.y) * v + rdUV.x * n;
}
int getSourceID()
{
	if (sourcesNum == 0)return -1;
	float nd = (sa = lightSources[sourcesNum - 1].a) *
		random(gl_FragCoord.yx * times / vec2(size));
	int l = 0, r = int(sourcesNum) - 1;
	if (nd < lightSources[l].a)return 0;
	while (l != r)
		if (nd < lightSources[(l + r) >> 1].a)r = (l + r) >> 1;
		else l = ((l + r) >> 1) + 1;
	return l;
}
float getPlaneT(Ray ray, vec4 plane)
{
	return -dot(plane, ray.p0) / dot(plane.xyz, ray.n);
}
vec2 getTriangleUV(vec3 pos, uint num)
{

	vec3 d = pos - triangles[num].p1;
	return vec2(dot(d, triangles[num].k1), dot(d, triangles[num].k2));
}
bool triangleTest(vec2 uv)
{
	return  all(greaterThanEqual(uv, vec2(0, 0))) && (uv.x + uv.y <= 1);
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
bool judgeHit(Ray ray)
{
	uint now = 0;
	uint bvhStack = 0;
	uint bvhState;
	int bvhSP = -1;
	for (;;)
	{
		if (judgeHitBox(ray, bvh[now].bound))
		{
			if (bvh[now].geometry != 0)
			{
				uint n = bvh[now].geometryNum;
				switch (bvh[now].geometry)
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
						if (tt > 0 && tt < ray.t)
						{
							vec2 uv = getTriangleUV(ray.p0.xyz + ray.n * tt, n);
							if (triangleTest(uv))return true;
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
							if (tt > 0 && tt < ray.t)return true;
						}
						break;
					}
					case 4:
					{
						float tt = getPlaneT(ray, circles[n].plane);
						if (tt > 0 && tt < ray.t)
						{
							vec3 d = ray.p0.xyz + ray.n * tt - circles[n].sphere.xyz;
							if (dot(d, d) <= circles[n].sphere.w)return true;
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
								if (tt > 0 && tt < ray.t)return true;
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
									tt = ttt;
							}
							if (tt > 0 && tt < ray.t)return true;
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
	return false;
}
HitInfo getFirstHit(Ray ray)
{
	uvec2 hitObj;
	ray.t = -1;
	uint now = 0;
	uint bvhStack = 0;
	uint bvhState;
	int bvhSP = -1;
	vec3 tempN;
	vec2 tempUV;
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
								hitObj = uvec2(geometry, n);
								//if (triangles[n].nIndices.x < 0)
								tempN = triangles[n].plane.xyz;
								//else
								//	tempN = normalize((1 - uv.x - uv.y) * water[triangles[n].nIndices.x].n +
								//		uv.x * water[triangles[n].nIndices.y].n +
								//		uv.y * water[triangles[n].nIndices.z].n);
								tempUV = (1 - uv.x - uv.y) * triangles[n].uv1 + uv.x * triangles[n].uv2 + uv.y * triangles[n].uv3;
							}
						}
						break;
					}
					case 3:
					{
						vec3 d = spheres[n].sphere.xyz - ray.p0.xyz;
						vec3 c = cross(d, ray.n);
						float s = spheres[n].sphere.w - dot(c, c);
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
								hitObj = uvec2(geometry, n);
								tempN = (ray.p0.xyz + ray.t * ray.n - spheres[n].sphere.xyz) / sqrt(spheres[n].sphere.w);
								float ne1 = dot(tempN, spheres[n].e1);
								float u =
									dot(tempN, cross(spheres[n].e1, spheres[n].e2)) >= 0 ?
									acos(dot(spheres[n].e2, tempN) / sqrt(1 - pow(ne1, 2))) / (2 * Pi) :
									1 - acos(dot(spheres[n].e2, tempN) / sqrt(1 - pow(ne1, 2))) / (2 * Pi);
								tempUV = vec2(u, 1 - acos(ne1) / Pi);
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
								hitObj = uvec2(geometry, n);
								tempN = circles[n].plane.xyz;
								vec3 e2 = cross(tempN, circles[n].e1);
								tempUV = (vec2(1) + vec2(dot(circles[n].e1, d), dot(e2, d)) / sqrt(circles[n].sphere.w)) / 2;
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
									hitObj = uvec2(geometry, n);
									tempN = normalize(d + ray.n * ray.t - cylinders[n].n * v);
									vec3 e2 = cross(cylinders[n].n, cylinders[n].e1);
									float u =
										dot(tempN, e2) >= 0 ?
										acos(dot(cylinders[n].e1, tempN)) / (2 * Pi) :
										1 - acos(dot(cylinders[n].e1, tempN)) / (2 * Pi);
									tempUV = vec2(u, v / cylinders[n].l);
								}
							}
						}
						break;
					}
					case 6:
					{
						vec3 d = ray.p0.xyz - cones[n].c;
						float nn0 = dot(ray.n, cones[n].n);
						float dn = dot(d, cones[n].n);
						vec3 j = dn * cones[n].n - cones[n].c2 * d;
						float c = dot(d, j);
						float b = dot(ray.n, j);
						float a = nn0 * nn0 - cones[n].c2;
						float s = b * b - a * c;
						if (s > 0)
						{
							s = sqrt(s);
							float tt = -1;
							float r2;
							tt = (s - b) / a;
							float d2 = dot(d, d);
							float dn0 = dot(d, ray.n);
							if (tt > 0)
							{
								r2 = d2 + tt * tt + 2 * dn0 * tt;
								float k = dn + nn0 * tt;
								if (r2 > cones[n].l2 || k < 0)tt = -1;
							}
							float ttt = (-b - s) / a;
							if (ttt > 0 && (ttt < tt || (tt < 0)))
							{
								float r2t = d2 + ttt * ttt + 2 * dn0 * ttt;
								float k = dn + nn0 * ttt;
								if (r2t <= cones[n].l2 && k > 0)
								{
									tt = ttt;
									r2 = r2t;
								}
							}
							if (tt > 0 && (tt < ray.t || ray.t < 0))
							{
								ray.t = tt;
								hitObj = uvec2(geometry, n);
								tempN = ((d + ray.n * ray.t) * sqrt(cones[n].c2 / r2) - cones[n].n) /
									sqrt(1 - cones[n].c2);
								vec3 nxy = normalize(d + ray.n * ray.t - cones[n].n * sqrt(r2 * cones[n].c2));
								float u =
									dot(nxy, cross(cones[n].n, cones[n].e1)) >= 0 ?
									acos(dot(cones[n].e1, nxy)) / (2 * Pi) :
									1 - acos(dot(cones[n].e1, nxy)) / (2 * Pi);
								tempUV = vec2(u, 1 - sqrt(r2 / cones[n].l2));
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
	return HitInfo(tempN, ray.t, tempUV, hitObj);
}
void forwardRayTrace()
{
	int sourceID = getSourceID();
	if (sourceID == -1)return;
	Ray ray;
	vec2 randUV = vec2(
		random((sqrt(gl_FragCoord.yx) + vec2(sqrt(times)))),
		random((sqrt(gl_FragCoord.xy) - vec2(times)))
	);
	uvec2 id = lightSources[sourceID].index;
	switch (id.x)
	{
		case 2:
		{
			float s = 1 - randUV.x - randUV.y;
			if (s < 0)
			{
				randUV = vec2(1) - randUV;
				s = -s;
			}
			subSources[0].color = trianglesOrigin[id.y].color.g * sa;
			float ahh = 1;
			while (subSources[0].color == vec3(0))
			{
				randUV = vec2(
					random((sqrt(gl_FragCoord.yx) + vec2(sqrt(times + ahh)))),
					random((sqrt(gl_FragCoord.xy) - vec2(times - ahh)))
				);
				s = 1 - randUV.x - randUV.y;
				if (s < 0)
				{
					randUV = vec2(1) - randUV;
					s = -s;
				}
				subSources[0].color = trianglesOrigin[id.y].color.g * sa;
				++ahh;
			}
			int tex = trianglesOrigin[id.y].color.texG;
			if (tex >= 0)
			{
				vec2 uv = s * trianglesOrigin[id.y].uv1 +
					randUV.x * trianglesOrigin[id.y].uv2 +
					randUV.y * trianglesOrigin[id.y].uv3;
				subSources[0].color *= texture(texSmp, vec3(uv, tex)).xyz;
			}
			subSources[0].n = (random(randUV) > 0.5 ? 1 : -1) * triangles[id.y].plane.xyz;
			ray.p0 = (subSources[0].p = vec4(
				s * trianglesOrigin[id.y].p1 +
				randUV.x * trianglesOrigin[id.y].p2 +
				randUV.y * trianglesOrigin[id.y].p3 +
				subSources[0].n * offset, 1));
			ray.n = getDirection(subSources[0].n,
				vec2(random(gl_FragCoord.yx + vec2(times)),
					random(gl_FragCoord.xy - vec2(sqrt(times)))));
			break;
		}
		case 3:
		case 4:
		case 5:
		case 6:
			break;
	}
	if (subSources[0].color == vec3(0))return;
	vec3 colorNow = subSources[0].color * dot(subSources[0].n, ray.n);
	ForwardStack forwardStack[ForwardDepth];
	int sp = -1;
	int depth = 0;
	HitInfo hitInfo;
	vec3 decayNow = decayOrigin.xyz;//to be changed
	subID = 1;
	while (subID < SubSoureMax)
	{
		hitInfo = getFirstHit(ray);
		Color tempColor;
		tempColor.texG = -1;
		if (hitInfo.t > 0 && depth < ForwardDepth)
		{
			switch (hitInfo.hitObj.x)
			{
				case 2:tempColor = triangles[hitInfo.hitObj.y].color; break;
				case 3:tempColor = spheres[hitInfo.hitObj.y].color; break;
				case 4:tempColor = circles[hitInfo.hitObj.y].color; break;
				case 5:tempColor = cylinders[hitInfo.hitObj.y].color; break;
				case 6:tempColor = cones[hitInfo.hitObj.y].color; break;
			}
			colorNow *= exp(decayNow * hitInfo.t);
			if (tempColor.texR >= 0)
				tempColor.r *= texture(texSmp, vec3(hitInfo.uv, tempColor.texR)).xyz;
			if (tempColor.texT >= 0)
				tempColor.t *= texture(texSmp, vec3(hitInfo.uv, tempColor.texT)).xyz;
			if (tempColor.texD >= 0)
				tempColor.d *= texture(texSmp, vec3(hitInfo.uv, tempColor.texD)).xyz;
			if (tempColor.texG >= 0)
				tempColor.g *= texture(texSmp, vec3(hitInfo.uv, tempColor.texG)).xyz;
			colorNow += tempColor.g;
			tempColor.t *= colorNow;
			float cosi1 = dot(ray.n, hitInfo.n);
			bool refracted = false;
			if (any(greaterThanEqual(tempColor.t, vec3(minColor))))
			{
				if (cosi1 > 0) tempColor.n = 1 / tempColor.n;
				float sini1 = sqrt(1 - cosi1 * cosi1);
				float sini2 = sini1 / tempColor.n;
				if (sini2 < 1)
				{
					float cosi2 = sqrt(1 - sini2 * sini2);
					if (sini2 <= 0.01)
					{
						float nadd1 = 1 / (tempColor.n + 1);
						tempColor.r *= pow((tempColor.n - 1) * nadd1, 2);
						tempColor.t *= pow(2 * nadd1, 2) * tempColor.n;
					}
					else
					{
						float a1 = tempColor.n * abs(cosi1) + cosi2;
						float a2 = abs(cosi1) + tempColor.n * cosi2;
						tempColor.r *= (pow((tempColor.n * cosi2 - abs(cosi1)) / a2, 2) +
							pow((cosi2 - tempColor.n * abs(cosi1)) / a1, 2)) / 2;
						tempColor.t *= 2 * cosi2 * (1 / pow(a1, 2) + 1 / pow(a2, 2)) * tempColor.n * abs(cosi1);
					}
					if (any(greaterThanEqual(tempColor.t, vec3(minColor))))
					{
						refracted = true;
						forwardStack[++sp].tDecayFactor = decayNow - sign(cosi1) * tempColor.decayFactor;
						forwardStack[sp].tp0 = ray.p0 + vec4((hitInfo.t + offset) * ray.n, 0);
						forwardStack[sp].tn = (ray.n + (tempColor.n * sign(cosi1) * cosi2 - cosi1) * hitInfo.n) / tempColor.n;
						forwardStack[sp].t = 1;
						forwardStack[sp].tColor = tempColor.t;
						forwardStack[sp].depth = depth + 1;
						forwardStack[sp].d = 0;
					}
				}
				else
				{
					tempColor.r = vec3(1);
				}
			}
			if (tempColor.d != vec3(0))
			{
				subSources[subID].p = ray.p0;
				subSources[subID].n = -sign(cosi1) * hitInfo.n;
				ray.p0 += vec4((hitInfo.t - offset) * ray.n, 0);
				uint n = 0;
				vec3 ts = colorNow * abs(cosi1);
				for (; n < pointLightNum; ++n)
				{
					vec3 dn = pointLights[n].p - ray.p0.xyz;
					float tt = dot(dn, dn);
					dn = normalize(dn);
					float ds = sign(cosi1) * dot(hitInfo.n, dn);
					if (ds < 0)
						if (!judgeHit(Ray(ray.p0, dn, sqrt(tt))))
						{
							ts -= (ds / tt) * pointLights[n].color;
						}
				}
				ts *= tempColor.d;
				subSources[subID++].color = ts;
				if (subID == SubSoureMax)return;
				vec3 u;
				if (!refracted)
				{
					forwardStack[++sp].depth = depth + 1;
					forwardStack[sp].t = 0;
				}
				forwardStack[sp].dp0 = ray.p0;
				forwardStack[sp].dn =
					getDirection(sign(cosi1) * hitInfo.n, vec2(
						random(vec2(sqrt(times)) + sqrt(gl_FragCoord.yx) + ray.n.xz * ray.n.yz),
						random(vec2(times) - sqrt(gl_FragCoord.xy) - ray.n.xz * ray.n.yx)));
				forwardStack[sp].d = 1;
				forwardStack[sp].dColor = ts;
				forwardStack[sp].dDecayFactor = decayNow;
			}
			colorNow *= tempColor.r;
			if (any(greaterThanEqual(colorNow, vec3(minColor))))
			{
				ray.n -= 2 * cosi1 * hitInfo.n;
				++depth;
				continue;
			}
		}
		if (sp < 0)return;
		else
		{
			depth = forwardStack[sp].depth;
			if (forwardStack[sp].t == 1)
			{
				ray.p0 = forwardStack[sp].tp0;
				ray.n = forwardStack[sp].tn;
				colorNow = forwardStack[sp].tColor;
				decayNow = forwardStack[sp].tDecayFactor;
				forwardStack[sp].t = 0;
				if (forwardStack[sp].d == 0)
					--sp;
				continue;
			}
			else if (forwardStack[sp].d == 1)
			{
				ray.p0 = forwardStack[sp].dp0;
				ray.n = forwardStack[sp].dn;
				colorNow = forwardStack[sp].dColor;
				decayNow = forwardStack[sp].dDecayFactor;
				forwardStack[sp].d = 0;
			}
			--sp;
		}
	}
	return;
}
Ray backwardRayAlloctor(float x)
{
	return Ray
	(
		vec4(r0, 1),
		normalize(trans * vec3(2 * (gl_FragCoord.xy +
			vec2(random(vec2(times + 6.275432 * x)), random(vec2(sqrt(times + 3.50975 * x))))
			) - size, z0)), -1
	);
}
vec4 backwardRayTrace(Ray ray)
{
	BackwardStack stack[BackwardDepth];
	int sp = -1;
	int depth = 0;
	vec3 ratioNow = vec3(1);
	vec3 answer = vec3(0);
	vec3 tempAnswer = vec3(0);
	HitInfo hitInfo;
	vec3 decayNow = decayOrigin.xyz;
	uint weight = 1;
	bool hittedSource = false;
	for (;;)
	{
		hitInfo = getFirstHit(ray);
		Color tempColor;
		tempColor.texG = -1;
		if (hitInfo.hitObj.x != 0)
		{
			switch (hitInfo.hitObj.x)
			{
				case 2:tempColor = triangles[hitInfo.hitObj.y].color; break;
				case 3:tempColor = spheres[hitInfo.hitObj.y].color; break;
				case 4:tempColor = circles[hitInfo.hitObj.y].color; break;
				case 5:tempColor = cylinders[hitInfo.hitObj.y].color; break;
				case 6:tempColor = cones[hitInfo.hitObj.y].color; break;
			}
		}
		if (tempColor.texG >= 0)
			tempColor.g *= texture(texSmp, vec3(hitInfo.uv, tempColor.texG)).xyz;
		if (hitInfo.t < 0)
		{
			hittedSource = true;
			tempColor.g = vec3(0);// texture(cubeSmp, ray.n).xyz;
			hitInfo.t = 0;
		}
		ratioNow *= exp(decayNow * hitInfo.t);
		if (tempColor.g != vec3(0))
			hittedSource = true;
		tempAnswer += tempColor.g * ratioNow;
		if (hitInfo.t > 0 && depth < BackwardDepth)
		{
			if (tempColor.texR >= 0)
				tempColor.r *= texture(texSmp, vec3(hitInfo.uv, tempColor.texR)).xyz;
			if (tempColor.texT >= 0)
				tempColor.t *= texture(texSmp, vec3(hitInfo.uv, tempColor.texT)).xyz;
			if (tempColor.texD >= 0)
				tempColor.d *= texture(texSmp, vec3(hitInfo.uv, tempColor.texD)).xyz;
			tempColor.t *= ratioNow;
			tempColor.d *= ratioNow;
			float cosi1 = dot(ray.n, hitInfo.n);
			bool refracted = false;
			if (any(greaterThanEqual(tempColor.t, vec3(minColor))))
			{
				if (cosi1 > 0) tempColor.n = 1 / tempColor.n;
				float sini1 = sqrt(1 - cosi1 * cosi1);
				float sini2 = sini1 / tempColor.n;
				if (sini2 < 1)
				{
					float cosi2 = sqrt(1 - sini2 * sini2);
					if (sini2 <= 0.01)
					{
						float nadd1 = 1 / (tempColor.n + 1);
						tempColor.r *= pow((tempColor.n - 1) * nadd1, 2);
						tempColor.t *= pow(2 * nadd1, 2) * tempColor.n;
					}
					else
					{
						float a1 = tempColor.n * abs(cosi1) + cosi2;
						float a2 = abs(cosi1) + tempColor.n * cosi2;
						tempColor.r *=
							(
								pow((cosi2 - tempColor.n * abs(cosi1)) / a1, 2) +
								pow((tempColor.n * cosi2 - abs(cosi1)) / a2, 2)
								) / 2;
						tempColor.t *= 2 * cosi2 * (1 / pow(a1, 2) + 1 / pow(a2, 2)) * tempColor.n * abs(cosi2);
					}
					if (any(greaterThanEqual(tempColor.t, vec3(minColor))))
					{
						refracted = true;
						stack[++sp].tDecayFactor = decayNow - sign(cosi1) * tempColor.decayFactor;
						stack[sp].tp0 = ray.p0 + vec4((hitInfo.t + offset) * ray.n, 0);
						stack[sp].tn = (ray.n + (tempColor.n * sign(cosi1) * cosi2 - cosi1) * hitInfo.n) / tempColor.n;
						stack[sp].t = 1;
						stack[sp].tRatio = tempColor.t;
						stack[sp].depth = depth + 1;
						stack[sp].d = 0;
					}
				}
				else
				{
					tempColor.r = vec3(1);
				}
			}
			ray.p0 += vec4((hitInfo.t - offset) * ray.n, 0);
			if (any(greaterThanEqual(tempColor.d, vec3(minColor))))
			{
				uint n;
				for (n = 0; n < pointLightNum; ++n)
				{
					vec3 dn = pointLights[n].p - ray.p0.xyz;
					float tt = dot(dn, dn);
					dn /= sqrt(tt);
					float ds = sign(cosi1) * dot(hitInfo.n, dn);
					if (ds < 0)
						if (!judgeHit(Ray(ray.p0, dn, sqrt(tt))))
							tempAnswer -= (ds / tt) * pointLights[n].color * tempColor.d;
				}
				vec3 ks = vec3(0);
				if (subID != 0)
				{
					vec3 dn;
					float ns;
					/*= subSources[0].p.xyz - ray.p0.xyz;
					float tt = dot(dn, dn);
					dn /= sqrt(tt); = dot(dn, subSources[0].n);
					if (ns < 0)
					{
						dn += subSources[0].n * offset;
						float ds = sign(cosi1) * dot(hitInfo.n, dn);
						if (ds < 0)
							if (!judgeHit(Ray(ray.p0, dn, sqrt(tt))))
							{
								ks += (ds * ns / tt) * subSources[0].color;
								++weight;
							}
					}*/
					for (n = 1; n < subID; ++n)
					{
						dn = subSources[n].p.xyz - ray.p0.xyz;
						float tt = dot(dn, dn);
						dn /= sqrt(tt);
						ns = dot(dn, normalize(subSources[n].n));
						if (ns < 0)
						{
							dn += subSources[n].n * offset;
							if (tt < 0.05)continue;
							float ds = sign(cosi1) * dot(hitInfo.n, dn);
							if (ds < 0)
								if (!judgeHit(Ray(ray.p0, dn, sqrt(tt))))
								{
									ks += (ds * ns) * subSources[n].color;
								}
						}
					}
					answer += ks * tempColor.d;
					if (ks != vec3(0))++weight;
				}
				vec3 u;
				float cos2Theta = random((gl_FragCoord.xy + vec2(sqrt(times)) * ray.p0.zx) / vec2(size));
				vec3 v;
				float phi = 2.0 * Pi * random((gl_FragCoord.xy + vec2(times) * ray.p0.yz) / vec2(size));
				if (abs(hitInfo.n.x) > 0.7)
				{
					float s = sqrt(1 - hitInfo.n.y * hitInfo.n.y);
					u = vec3(-hitInfo.n.z, 0, hitInfo.n.x) / s;
					v = vec3(hitInfo.n.x * hitInfo.n.y / s, -s,
						hitInfo.n.y * hitInfo.n.z / s);
				}
				else
				{
					float s = sqrt(1 - hitInfo.n.x * hitInfo.n.x);
					u = vec3(0, hitInfo.n.z, -hitInfo.n.y) / s;
					v = vec3(-s, hitInfo.n.x * hitInfo.n.y / s,
						hitInfo.n.x * hitInfo.n.z / s);
				}
				if (!refracted)
				{
					stack[++sp].depth = depth + 1;
					stack[sp].t = 0;
				}
				float sinTheta = sqrt(1 - cos2Theta);
				float cc = sqrt(cos2Theta);
				stack[sp].dp0 = ray.p0;
				stack[sp].dn =
					sinTheta * cos(phi) * u +
					sinTheta * sin(phi) * v -
					sign(cosi1) * cc * hitInfo.n;
				stack[sp].d = 1;
				stack[sp].dRatio = tempColor.d;
				stack[sp].dDecayFactor = decayNow;
			}
			ratioNow *= tempColor.r;
			if (any(greaterThanEqual(ratioNow, vec3(minColor))))
			{
				ray.n -= 2 * cosi1 * hitInfo.n;
				++depth;
				continue;
			}
		}
		if (sp < 0)
		{
			if (tempAnswer == vec3(0) && !hittedSource)--weight;
			else answer += tempAnswer;
			return vec4(answer, weight);
		}
		else
		{
			depth = stack[sp].depth;
			if (stack[sp].t == 1)
			{
				ray.p0 = stack[sp].tp0;
				ray.n = stack[sp].tn;
				ratioNow = stack[sp].tRatio;
				decayNow = stack[sp].tDecayFactor;
				stack[sp].t = 0;
				if (stack[sp].d == 0)
					--sp;
				continue;
			}
			else if (stack[sp].d == 1)
			{
				ray.p0 = stack[sp].dp0;
				ray.n = stack[sp].dn;
				ratioNow = stack[sp].dRatio;
				decayNow = stack[sp].dDecayFactor;
				stack[sp].d = 0;
			}
			--sp;
		}
	}
	if (tempAnswer == vec3(0) && !hittedSource)--weight;
	else answer += tempAnswer;
	return vec4(answer, weight);
}

void main()
{
	subID = 0;
	forwardRayTrace();
	vec4 rt = backwardRayTrace(backwardRayAlloctor(1));
	rt += backwardRayTrace(backwardRayAlloctor(2));
	//atomicAdd(counter, 1);
	vec4 c = imageLoad(image, ivec2(gl_FragCoord.xy));
	if (times == 1)
	{
		imageStore(image, ivec2(gl_FragCoord.xy), rt);
		if (rt.w == 0)
			gl_FragColor = rt;
		else
			gl_FragColor = vec4(rt.xyz / rt.w, rt.w);
	}
	else
	{
		imageStore(image, ivec2(gl_FragCoord.xy), rt += c);
		if (c.w + rt.w == 0)
			gl_FragColor = vec4(0);
		else
			gl_FragColor = vec4(rt.xyz / rt.w, rt.w);
	}
}
