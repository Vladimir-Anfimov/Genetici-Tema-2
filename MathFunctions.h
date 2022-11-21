#pragma once
#include <vector>
#include <corecrt_math_defines.h>

using std::vector;

enum MATH_FUNC {
	DE_JONG = 1,
	MICHALEWICS = 2,
	RASTRIGIN = 3,
	SCHWEFEL = 4
};

class MathFunctions
{
public:
	static float de_jong(const vector<float>& x_real)
	{
		float sum = 0;
		for (int i = 0; i < x_real.size(); i++)
		{
			sum += pow(x_real[i], 2);
		}
		return sum;
	}
	static float Michalewics(const vector<float>& x_real)
	{
		float sum = 0;
		for (int i = 0; i < x_real.size(); i++)
		{
			sum += sin(x_real[i]) * pow(sin((i + 1) * x_real[i] * x_real[i] / M_PI), 20);
		}
		return (-1) * sum;
	}

	static float Rastrigin(const vector<float>& x_real)
	{
		float sum = 10 * x_real.size();
		for (int i = 0; i < x_real.size(); i++)
		{
			sum += x_real[i] * x_real[i] - 10 * cos(2 * M_PI * x_real[i]);
		}
		return sum;
	}

	static double Schwefel(const vector<float>& x_real)
	{
		float sum = 0;
		for (int i = 0; i < x_real.size(); i++)
		{
			sum -= x_real[i] * sin(sqrt(abs(x_real[i])));
		}
		return sum;
	}
};


