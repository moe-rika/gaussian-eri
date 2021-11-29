#include "test_by_simpson.h"

#include <iostream>
#include<array>
#include<cmath>

using std::array;
using std::cout;
using std::endl;

#define Power pow

template <unsigned int N>
struct HighDimSimpsonIntegration
{
	struct IntegrationDescriptor
	{
		int div_;
		double min, max, val;
	};
	static double simpson_int(
		std::array<IntegrationDescriptor, N>& d,
		unsigned int i,
		double(*f)(std::array<IntegrationDescriptor, N>&)
	)
	{

		int& div_ = d[i].div_;
		double& min = d[i].min, & max = d[i].max, & val = d[i].val;
		if (i != N - 1)
		{
			double sum = 0;
			val = min;
			sum += simpson_int(d, i + 1, f);
			val = max;
			sum += simpson_int(d, i + 1, f);
			double step = (max - min) / div_;
			for (int j = 1; j < div_; j += 2)
			{
				val = min + step * j;
				sum += 4 * simpson_int(d, i + 1, f);
				val = min + step * (j + 1);
				sum += 2 * simpson_int(d, i + 1, f);
				if (i == 0)
					std::cout << j << std::endl;
			}
			return (max - min) / 3 / div_ * sum;
		}
		else
		{
			double sum = 0;
			val = min;
			sum += f(d);
			val = max;
			sum += f(d);
			double step = (max - min) / div_;
			for (int j = 1; j < div_; j += 2)
			{
				val = min + step * j;
				sum += 4 * f(d);
				val = min + step * (j + 1);
				sum += 2 * f(d);
			}
			return (max - min) / 3 / div_ * sum;
		}
	}
};

double f(double c, double a,
	int ax, int ay, int az,
	double x0, double y0, double z0,
	double x, double y, double z)
{
	return (c * Power(x - x0, ax) * Power(y - y0, ay) * Power(z - z0, az)) / exp(a * (Power(x - x0, 2) + Power(y - y0, 2) + Power(z - z0, 2)));
}

inline double f1(double x, double y, double z)
{
	return f(1, 0.12, 1, 0, 0, 0, 0, 0, x, y, z);
}

inline double f2(double x, double y, double z)
{
	return f(1, 0.23, 1, 0, 0, 1, 0, 0, x, y, z);
}

inline double f3(double x, double y, double z)
{
	return f(1, 0.34, 0, 2, 0, 0, 0, 0, x, y, z);
}

inline double f4(double x, double y, double z)
{
	return f(1, 0.45, 0, 2, 0, 1, 0, 0, x, y, z);
}

double func(std::array<HighDimSimpsonIntegration<6>::IntegrationDescriptor, 6>& d)
{
	const auto& x1 = d[0].val;
	const auto& y1 = d[1].val;
	const auto& z1 = d[2].val;
	const auto& x2 = d[3].val;
	const auto& y2 = d[4].val;
	const auto& z2 = d[5].val;
	return (f1(x1, y1, z1) * f2(x1, y1, z1) * f3(x2, y2, z2) * f4(x2, y2, z2)) /
		sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
}




const int div_ = 25, dim = 6;
const double range = 6;
std::array<HighDimSimpsonIntegration<dim>::IntegrationDescriptor, dim> d{
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range + range / div_ / 2,range + range / div_ / 2,0},
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range + range / div_ / 2,range + range / div_ / 2,0},
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range + range / div_ / 2,range + range / div_ / 2,0},
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range,range,0},
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range,range,0},
	HighDimSimpsonIntegration<dim>::IntegrationDescriptor
{div_,-range,range,0}
};


int simpson_main()
{
	double k = HighDimSimpsonIntegration<6>::simpson_int(d, 0, func);
	std::cout << k << std::endl;//81.7827
	return 0;
}
