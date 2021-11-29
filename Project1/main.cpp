#include<vector>
#include<deque>
#include<iostream>
#include<array>
#include<cmath>

#include "Point3D.hpp"

using std::vector;
using std::deque;
using std::cout;
using std::endl;
using std::array;

using AngularMomentum = array<signed char, 3>;


const double PI = 3.141592653589793238463;

class BasisSet
{
public:
	AngularMomentum a;
	double coeficient;
	double alpha;
	Point3D center;
};

class BasisSetFourTuple
{
public:
	BasisSetFourTuple(const BasisSet& a, const BasisSet& b, const BasisSet& c, const BasisSet& d) :
		coeff{ a.coeficient,b.coeficient,c.coeficient,d.coeficient }
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
		a_ = a.a;
		b_ = b.a;
		c_ = c.a;
		d_ = d.a;
		P = (alpha * A + beta * B) / xi;
		Q = (gamma * C + delta * D) / eta;
		T = rho * (P - Q).norm_sq();
		Sab = exp(-alpha * beta / xi * (A - B).norm_sq());
		Scd = exp(-gamma * delta / eta * (C - D).norm_sq());

		A_B = A - B;
		C_D = C - D;
		RR = Q - C + xi / eta * (P - A);
		two_eta = eta * 2;
		neg_xi_div_eta = -xi / eta;
		P_A = P - A;
		P_Q = -rho / xi * (P - Q);
		two_xi = 2 * xi;
		neg_rho_div_xi = -rho / xi;
		coef = pow(3.141592657 / (xi + eta), 1.5) * Sab * Scd;
		dd.push_back(BasisSetFourTupleUnit{ a_, b_, c_, d_, 1 });
	}


	vector<double> coeff;
	AngularMomentum a_, b_, c_, d_;
	double alpha, beta, gamma, delta;
	double xi, eta, rho, T, Sab, Scd;
	Point3D A, B, C, D, P, Q;


	double two_eta, neg_xi_div_eta, two_xi, neg_rho_div_xi;
	Point3D A_B, C_D, RR, P_A, P_Q;
	double coef;


	struct BasisSetFourTupleUnit
	{
		AngularMomentum a, b, c, d;
		double coeff;
	};

	struct BasisSetFourTupleUnit1
	{
		AngularMomentum a, c, d;
		double coeff;
	};

	struct BasisSetFourTupleUnit2
	{
		AngularMomentum a, c;
		double coeff;
	};

	struct BasisSetFourTupleUnit3
	{
		AngularMomentum a;
		signed char m;
		double coeff;
	};

	struct BasisSetFourTupleUnit4
	{
		signed char m;
		double coeff;
	};

	deque<BasisSetFourTupleUnit> dd;
	deque<BasisSetFourTupleUnit1> dd1;
	deque<BasisSetFourTupleUnit2> dd2;
	deque<BasisSetFourTupleUnit3> dd3;
	deque<BasisSetFourTupleUnit4> dd4;



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


};



int main()
{
	BasisSet px1{ {1,0,0},1,0.12,{0,0,0} };
	BasisSet px2{ {1,0,0},1,0.23,{1,0,0} };
	BasisSet px3{ {0,2,0},1,0.34,{0,0,0} };
	BasisSet px4{ {0,2,0},1,0.45,{1,0,0} };
	auto m0 = (BasisSetFourTuple(px1, px2, px3, px4));
	auto m = m0;
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


void BasisSetFourTuple::HRR1() //(ab|cd) -> (a0|cd)
{
	deque<BasisSetFourTupleUnit> temp;
	for (signed char i = 0; i < 3; i++)
	{
		while (!dd.empty())
		{
			auto& front = dd.front();
			if (front.b[i] == 0)
			{
				temp.push_back(front);
				dd.pop_front();
			}
			else
			{

			}
		}
	}
	//for (int i = 0; i < 3; i++)
	//{
	//	auto iter = m.begin();
	//	while (iter != m.end()) {
	//		if (iter->first[3 + i] == 0)
	//		{
	//			++iter;
	//		}
	//		else
	//		{
	//			vector<signed char> first = iter->first;
	//			vector<signed char> temp = iter->first;
	//			double second = iter->second;
	//			m.erase(iter);

	//			temp[i] = first[i] + 1;
	//			temp[3 + i] = first[3 + i] - 1;
	//			if (m.find(temp) == m.end())
	//			{
	//				m[temp] = second * 1.0;
	//			}
	//			else
	//			{
	//				m[temp] += second * 1.0;
	//			}

	//			temp = first;
	//			temp[i] = first[i];
	//			temp[3 + i] = first[3 + i] - 1;
	//			if (m.find(temp) == m.end())
	//			{
	//				m[temp] = second * A_B[i];
	//			}
	//			else
	//			{
	//				m[temp] += second * A_B[i];
	//			}
	//			iter = m.begin();
	//		}
	//	}
	//}
}
void BasisSetFourTuple::HRR2() //(a0|cd) -> (a0|c0)
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
void BasisSetFourTuple::HRR3() //(a0|c0) -> (a0|00)
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
void BasisSetFourTuple::VRR()//(a0|00) -> (00|00)
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


//0 0 0 0  1.446
//0 0 0 1 - 2.10049
//0 0 0 2  0.885293
//0 0 0 3 - 0.117013
//0 0 0 4  0.000416695
//0 0 0 5  0
//0 0 0 6  0
//234.811
//0 0.999381 0.00185794 1.44511
//1 0.332962 0.00185794 - 0.699382
//2 0.199735 0.00185794 0.176824
//3 0.142651 0.00185794 - 0.016692
//4 0.110892 0.00185794 4.62082e-05
//5 - 0.0308185 0.00185794 - 0
//6 - 359.848 0.00185794 - 0
//81.7487
