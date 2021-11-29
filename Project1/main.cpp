//Obara Daika shceme 
//DOI: 10.1039/b605188j

#include<vector>
#include<map>
#include<iostream>
#include<array>
#include<cmath>

using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::array;

const double PI = 3.141592653589793238463;

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

class BasisSet
{
public:
	signed char a[3];
	double coeficient;
	double alpha;
	Point3D center;
};

class BasisSetFourToulpe
{
public:
	BasisSetFourToulpe(const BasisSet& a, const BasisSet& b, const BasisSet& c, const BasisSet& d) :
		v{ a.a[0],a.a[1],a.a[2],b.a[0],b.a[1],b.a[2],c.a[0],c.a[1],c.a[2],d.a[0],d.a[1],d.a[2] },
		c{ a.coeficient,b.coeficient,c.coeficient,d.coeficient }
	{
		xi = a.alpha + b.alpha;
		eta = c.alpha + d.alpha;
		rho = xi * eta / (xi + eta);
		alpha = a.alpha;
		beta = b.alpha;
		gamma = c.alpha;
		delta = d.alpha;
		A = a.center;
		B = b.center;
		C = c.center;
		D = d.center;
		P = (alpha * A + beta * B) / xi;
		Q = (gamma * C + delta * D) / eta;
		T = rho * (P - Q).norm_sq();
		Sab = exp(-alpha * beta / xi * (A - B).norm_sq());
		Scd = exp(-gamma * delta / eta * (C - D).norm_sq());
	}
	vector<signed char> v;
	vector<double> c;
	double alpha, beta, gamma, delta;
	double xi, eta, rho, T, Sab, Scd;
	Point3D A, B, C, D, P, Q;
};

class BasisSetFourToulpeRecursionInterace
{
public:
	BasisSetFourToulpeRecursionInterace(const BasisSetFourToulpe& p)
	{
		m[p.v] = 1.0;
		A_B = p.A - p.B;
		C_D = p.C - p.D;
		RR = p.Q - p.C + p.xi / p.eta * (p.P - p.A);
		two_eta = p.eta * 2;
		neg_xi_div_eta = -p.xi / p.eta;
		P_A = p.P - p.A;
		P_Q = -p.rho / p.xi * (p.P - p.Q);
		two_xi = 2 * p.xi;
		neg_rho_div_xi = -p.rho / p.xi;
		coef = pow(3.141592657 / (p.xi + p.eta), 1.5) * p.Sab * p.Scd ;
	}




	void HRR1();//I(ab|cd) = I((a+1)(b-1)|cd) + (A_i - B_i)I(a(b-1)|cd)
	void HRR2();//I(ab|cd) = I(ab|(c+1)(d-1)) + (C_i - D_i)I(ab|c(d-1))
	void HRR3();
	void VRR();

	double boys_function(int n, double T)
	{
		if (n == 0)
		{
			return sqrt(PI / (4 * T)) * erf(sqrt(T));
		}
		else
		{
			return 1. / (2 * T) * ((2 * n - 1) * boys_function(n - 1, T) - exp(-T));
		}
	}

	map<vector<signed char>, double> m;
	map<vector<signed char>, double> m1;
	double two_eta, neg_xi_div_eta, two_xi, neg_rho_div_xi;
	Point3D A_B, C_D, RR, P_A, P_Q;
	double coef;
};



int main()
{
	BasisSet px1{ {1,0,0},1,0.12,{0,0,0} };
	BasisSet px2{ {1,0,0},1,0.23,{1,0,0} };
	BasisSet px3{ {0,2,0},1,0.34,{0,0,0} };
	BasisSet px4{ {0,2,0},1,0.45,{1,0,0} };
	auto m0 = (BasisSetFourToulpe(px1, px2, px3, px4));
	auto m = BasisSetFourToulpeRecursionInterace(m0);
	m.HRR1();
	m.HRR2();
	m.HRR3();
	m.VRR();
	for (auto& a : m.m1)
	{
		for (auto& b : a.first)
		{
			cout << (int)b << " ";
		}
		cout << " ";
		cout << a.second << endl;
	}
	double coeff;
	coeff = m.coef * m0.c[0] * m0.c[1] * m0.c[2] * m0.c[3];

	cout << coeff * pow(3.141592657 / m0.rho, 1.5) * m.m1.begin()->second << endl;
	double result = 0;
	for (auto& a : m.m1)
	{
		result += a.second * m.boys_function((int)a.first[3], m0.T);
		cout << (int)a.first[3] << " " << m.boys_function((int)a.first[3], m0.T) << " "
			<< m0.T << " "
			<< a.second * m.boys_function((int)a.first[3], m0.T)
			<< endl;
	}

	result *= coeff * 2 * PI / m0.rho;

	cout << result << endl;

	return 0;
}


void BasisSetFourToulpeRecursionInterace::HRR1() //(ab|cd) -> (a0|cd)
{
	for (int i = 0; i < 3; i++)
	{
		auto iter = m.begin();
		while (iter != m.end()) {
			if (iter->first[3 + i] == 0)
			{
				++iter;
			}
			else
			{
				vector<signed char> first = iter->first;
				vector<signed char> temp = iter->first;
				double second = iter->second;
				m.erase(iter);

				temp[i] = first[i] + 1;
				temp[3 + i] = first[3 + i] - 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * 1.0;
				}
				else
				{
					m[temp] += second * 1.0;
				}

				temp = first;
				temp[i] = first[i];
				temp[3 + i] = first[3 + i] - 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * A_B[i];
				}
				else
				{
					m[temp] += second * A_B[i];
				}
				iter = m.begin();
			}
		}
	}
}
void BasisSetFourToulpeRecursionInterace::HRR2() //(a0|cd) -> (a0|c0)
{
	for (int i = 0; i < 3; i++)
	{
		auto iter = m.begin();
		while (iter != m.end()) {
			if (iter->first[9 + i] == 0)
			{
				++iter;
			}
			else
			{
				vector<signed char> first = iter->first;
				vector<signed char> temp = iter->first;
				double second = iter->second;
				m.erase(iter);

				temp[6 + i] = first[6 + i] + 1;
				temp[9 + i] = first[9 + i] - 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * 1.0;
				}
				else
				{
					m[temp] += second * 1.0;
				}

				temp = first;
				temp[6 + i] = first[6 + i];
				temp[9 + i] = first[9 + i] - 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * C_D[i];
				}
				else
				{
					m[temp] += second * C_D[i];
				}
				iter = m.begin();
			}
		}
	}
}
void BasisSetFourToulpeRecursionInterace::HRR3() //(a0|c0) -> (a0|00)
{
	for (int i = 0; i < 3; i++)
	{
		auto iter = m.begin();
		while (iter != m.end()) {
			if (iter->first[6 + i] == 0)
			{
				++iter;
			}
			else
			{
				vector<signed char> first = iter->first;
				vector<signed char> temp = iter->first;
				double second = iter->second;
				m.erase(iter);

				temp[6 + i] = first[6 + i] - 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * RR[i];
				}
				else
				{
					m[temp] += second * RR[i];
				}

				if (first[i] != 0)
				{
					temp = first;
					temp[6 + i] = first[6 + i] - 1;
					temp[i] = first[i] - 1;
					if (m.find(temp) == m.end())
					{
						m[temp] = second * first[i] / two_eta;
					}
					else
					{
						m[temp] += second * first[i] / two_eta;
					}
				}

				if (first[6 + i] != 1)
				{
					temp = first;
					temp[6 + i] = first[6 + i] - 2;
					if (m.find(temp) == m.end())
					{
						m[temp] = second * (first[6 + i] - 1) / two_eta;
					}
					else
					{
						m[temp] += second * (first[6 + i] - 1) / two_eta;
					}
				}

				temp = first;
				temp[6 + i] = first[6 + i] - 1;
				temp[i] = first[i] + 1;
				if (m.find(temp) == m.end())
				{
					m[temp] = second * neg_xi_div_eta;
				}
				else
				{
					m[temp] += second * neg_xi_div_eta;
				}
				iter = m.begin();
			}
		}
	}
}
void BasisSetFourToulpeRecursionInterace::VRR()//(a0|00) -> (00|00)
{
	vector<signed char> v1;
	v1.resize(4); // a1 a2 a3 m
	for (auto& a : m)
	{
		for (int i = 0; i < 3; i++)
		{
			v1[i] = a.first[i];
		}
		v1[3] = 0;
		m1[v1] = a.second;
	}

	for (int i = 0; i < 3; i++)
	{
		auto iter = m1.begin();
		while (iter != m1.end()) {
			if (iter->first[i] == 0)
			{
				++iter;
			}
			else
			{
				vector<signed char> first = iter->first;
				vector<signed char> temp = iter->first;
				double second = iter->second;
				m1.erase(iter);

				temp[i] = first[i] - 1;
				if (m1.find(temp) == m1.end())
				{
					m1[temp] = second * P_A[i];
				}
				else
				{
					m1[temp] += second * P_A[i];
				}

				temp = first;
				temp[i] = first[i] - 1;
				temp[3] = first[3] + 1;
				if (m1.find(temp) == m1.end())
				{
					m1[temp] = second * P_Q[i];
				}
				else
				{
					m1[temp] += second * P_Q[i];
				}


				if (first[i] != 1)
				{
					temp = first;
					temp[i] = first[i] - 2;
					if (m1.find(temp) == m1.end())
					{
						m1[temp] = second * (first[i] - 1) / two_xi;
					}
					else
					{
						m1[temp] += second * (first[i] - 1) / two_xi;
					}

					temp = first;
					temp[i] = first[i] - 2;
					temp[3] = first[3] + 1;
					if (m1.find(temp) == m1.end())
					{
						m1[temp] = second * (first[i] - 1) / two_xi * neg_rho_div_xi;
					}
					else
					{
						m1[temp] += second * (first[i] - 1) / two_xi * neg_rho_div_xi;
					}
				}
				iter = m1.begin();
			}
		}
	}
}


