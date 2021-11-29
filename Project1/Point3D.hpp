#pragma once
struct Point3D
{
	double a[3];
	double& operator[](const int& b)
	{
		return a[b];
	}
	const double& operator[](const int& b) const
	{
		return a[b];
	}
	Point3D operator+(const Point3D& m)
	{
		return Point3D{ a[0] + m[0],a[1] + m[1],a[2] + m[2] };
	}
	Point3D operator-(const Point3D& m)
	{
		return Point3D{ a[0] - m[0],a[1] - m[1],a[2] - m[2] };
	}

	Point3D operator+(const Point3D& m) const
	{
		return Point3D{ a[0] + m[0],a[1] + m[1],a[2] + m[2] };
	}
	Point3D operator-(const Point3D& m) const
	{
		return Point3D{ a[0] - m[0],a[1] - m[1],a[2] - m[2] };
	}

	Point3D operator*(const double& m)
	{
		return Point3D{ a[0] * m,a[1] * m,a[2] * m };
	}

	friend Point3D operator*(const double& p, const Point3D& q)
	{
		return Point3D{ q[0] * p,q[1] * p,q[2] * p };
	}

	Point3D operator/(const double& m)
	{
		return Point3D{ a[0] / m,a[1] / m,a[2] / m };
	}

	double norm_sq()
	{
		return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
	}
};