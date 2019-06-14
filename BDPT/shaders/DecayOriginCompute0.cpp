#version 450 core
layout(local_size_x = 1)in;
#define originSamples 8


layout(std430, binding = 8)buffer DecayOrigin
{
	vec4 decayOrigins[originSamples];
	vec4 decayOrigin;
};

void main()
{
	bool judged[originSamples];
	int num[originSamples];
	int c0 = 0;
	for (; c0 < originSamples; ++c0)
	{
		judged[c0] = false;
		num[c0] = 0;
	}
	for (c0 = 0; c0 < originSamples; ++c0)
	{
		if (judged[c0])continue;
		judged[c0] = true;
		num[c0]++;
		int c1 = 0;
		for (; c1 < originSamples; ++c1)
		{
			if (!judged[c1] && decayOrigins[c1] == decayOrigins[c0])
			{
				num[c0]++;
				judged[c1] = true;
			}
		}
	}
	int index = 0;
	int max = 0;
	for (c0 = 0; c0 < originSamples; ++c0)
	{
		if (num[c0] > max)
		{
			max = num[c0];
			index = c0;
		}
	}
	decayOrigin = decayOrigins[index];
}