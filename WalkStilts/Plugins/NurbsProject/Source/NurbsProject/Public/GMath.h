// Copyright Zhangci 2018
//class of 3d math - vector 3d

#pragma once
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdarg>
#include <iomanip>
#include <vector>

using namespace std;

#define MAX_DEG 10

// pi
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923

// macro functions
#define SQRT(X)		sqrt((X))
#define SQR(X)		((X) * (X))
#define DEG2RAD(X)	((X) * (M_PI) / (180.0))
#define RAD2DEG(X)	((X) * (180.0) / (M_PI))
#define SWAP(type, x, y) { type temp = (x); (x) = (y); (y) = temp; }
#define MIN(x, y)	((x) > (y) ? (y) : (x))
#define MAX(x, y)	((x) > (y) ? (x) : (y))
#define ABS(X)		(((X) > 0.0) ? (X) : (-(X)))
#define SIGN(a)		((a) > 0.0 ? (1.0) : (-1.0))
#define SIGN1(a, b) ((b) > 0.0 ? ABS(a) : -ABS(a))
#define SIGN2(a, b)	((b) >= 0.0 ? fabs(a) : -fabs(a))
#define PYTHAG(a, b) SQRT((SQR(a) + SQR(b)))
#define EQ(X, Y, EPS)	(ABS((X) - (Y)) < EPS)
#define EQ_ZERO(X, EPS) (ABS(X) < EPS)
#define ARR_ZERO(A, N) memset((A), 0, sizeof(A[0]) * (N))
#define ARR_COPY(D, S, N) memmove((D), (S), sizeof(S[0]) * (N))
#define ERROR 0.000001


class GVector3
{
public:
	// Constructors and destructor
	GVector3(double x = 0.0, double y = 0.0, double z = 0.0);
	GVector3(const GVector3 &copy);
	virtual ~GVector3();

	// Assignment operator overloading
	GVector3 &operator =(const GVector3 &rhs);

	// Compound assignment operator overloading
	GVector3 &operator +=(const GVector3 &rhs);
	GVector3 &operator -=(const GVector3 &rhs);
	GVector3 &operator *=(const double &s);
	GVector3 &operator /=(const double &s);
	GVector3 &operator ^=(const GVector3 &rhs);

	// Unary operator overloading
	GVector3 operator +() const;
	GVector3 operator -() const;

	// Arithmetic operator overloading
	GVector3 operator +(const GVector3 &rhs) const;
	GVector3 operator -(const GVector3 &rhs) const;
	double   operator *(const GVector3 &rhs) const;
	GVector3 operator /(const double &s)	 const;
	GVector3 operator ^(const GVector3 &rhs) const;

	// Equality and inequality operator overloading
	bool operator ==(const GVector3 &rhs) const;
	bool operator !=(const GVector3 &rhs) const;

	// Subscript operator overloading
	double &operator [](const int &idx);
	const double &operator [](const int &idx) const;

	// Member functions
	GVector3 &Set(const double &x, const double &y, const double &z);
	GVector3 &Normalize();

	// Static member functions
	static void SetPrecision(double error);
	static double GetPrecision();

	// Declaration of friend
	friend GVector3 operator *(const GVector3 &lhs, const double &s);
	friend GVector3 operator *(const double &s, const GVector3 &rhs);
	friend ostream &operator <<(ostream &os, const GVector3 &v);
	friend istream &operator >> (istream &is, GVector3 &v);
	friend GVector3 proj(const GVector3 &v, const GVector3 &w);
	friend double dist(const GVector3 &v, const GVector3 &w);
	friend double norm(const GVector3 &v);
	friend double angle(const GVector3 &v, const GVector3 &w, bool radian);
	

protected:
	// Data members
	double V[3];
	static double Precision;
};

