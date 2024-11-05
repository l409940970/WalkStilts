// Copyright Zhangci 2018

#include "GMath.h"
#include "NurbsProject.h"


double GVector3::Precision = 0.00000001;


/*!
*	Default constructor
*
*	\param x x-coordinate
*	\param y y-coordinate
*	\param z z-coordinate
*/
GVector3::GVector3(double x, double y, double z)
{
	V[0] = x;
	V[1] = y;
	V[2] = z;
}

/*!
*	Copy constructor
*
*	\param copy GVector3 object to be copied
*/
GVector3::GVector3(const GVector3 &copy)
{
	V[0] = copy.V[0];
	V[1] = copy.V[1];
	V[2] = copy.V[2];
}

/*!
*	Destructor
*/
GVector3::~GVector3()
{
}

GVector3 &GVector3::operator =(const GVector3 &rhs)
{
	V[0] = rhs.V[0];
	V[1] = rhs.V[1];
	V[2] = rhs.V[2];
	return *this;
}

/*!
*	Compound assignment operator
*
*	\param rhs GVector3 object to be added
*/
GVector3 &GVector3::operator +=(const GVector3 &rhs)
{
	V[0] += rhs.V[0];
	V[1] += rhs.V[1];
	V[2] += rhs.V[2];
	return *this;
}

/*!
*	Compound assignment operator
*
*	\param rhs GVector3 object to be subtracted
*/
GVector3 &GVector3::operator -=(const GVector3 &rhs)
{
	V[0] -= rhs.V[0];
	V[1] -= rhs.V[1];
	V[2] -= rhs.V[2];
	return *this;
}

/*!
*	Compound assignment operator
*
*	\param s scaling factor
*/
GVector3 &GVector3::operator *=(const double &s)
{
	V[0] *= s;
	V[1] *= s;
	V[2] *= s;
	return *this;
}

/*!
*	Compound assignment operator
*
*	\param s division factor
*/
GVector3 &GVector3::operator /=(const double &s)
{
	V[0] /= s;
	V[1] /= s;
	V[2] /= s;
	return *this;
}

/*!
*	Compound assignment operator
*
*	\param rhs GVector3 object to be cross-producted
*/
GVector3 &GVector3::operator ^=(const GVector3 &rhs)
{
	double x = V[0], y = V[1], z = V[2];
	V[0] = y * rhs.V[2] - z * rhs.V[1];
	V[1] = z * rhs.V[0] - x * rhs.V[2];
	V[2] = x * rhs.V[1] - y * rhs.V[0];
	return *this;
}

/*!
*	Unary operator (+)
*/
GVector3 GVector3::operator +() const
{
	return *this;
}

/*!
*	Unary operator (-)
*/
GVector3 GVector3::operator -() const
{
	return *this * -1;
}

/*!
*	Operator + overloading
*
*	\param rhs right hand operand
*	\return a new vector that adds two vectors
*/
GVector3 GVector3::operator +(const GVector3 &rhs) const
{
	GVector3 ret(*this);
	ret += rhs;
	return ret;
}

/*!
*	Operator - overloading
*
*	\param rhs right hand operand
*	\return a new vector that subtracts two vectors
*/
GVector3 GVector3::operator -(const GVector3 &rhs) const
{
	GVector3 ret(*this);
	ret -= rhs;
	return ret;
}

/*!
*	Operator * (inner product)overloading
*
*	\param rhs right hand operand
*	\return the result of the inner product of two vectors
*/
double GVector3::operator *(const GVector3 &rhs) const
{
	return V[0] * rhs.V[0] + V[1] * rhs.V[1] + V[2] * rhs.V[2];
}

/*!
*	Operator / (divided)overloading
*
*	\param s scale factor
*	\return the vector divided by \a s
*/
GVector3 GVector3::operator /(const double &s) const
{
	GVector3 ret(*this);
	ret /= s;
	return ret;
}

/*!
*	Operator ^ (cross-produce) overloading
*
*	\param rhs right hand operand
*	\return the resulting vector of the cross product
*/
GVector3 GVector3::operator ^(const GVector3 &rhs) const
{
	return GVector3(V[1] * rhs.V[2] - V[2] * rhs.V[1], V[2] * rhs.V[0] - V[0] * rhs.V[2], V[0] * rhs.V[1] - V[1] * rhs.V[0]);
}

/*!
*	Equality operator overloading
*
*	\param rhs right hand operand
*	\return the result of equality test
*/
bool GVector3::operator ==(const GVector3 &rhs) const
{
	return !((*this) != rhs);
}

/*!
*	Inequality operator overloading
*
*	\param rhs right hand operand
*	\return the result of equality test
*/
bool GVector3::operator !=(const GVector3 &rhs) const
{
	return (!EQ(V[0], rhs.V[0], Precision) || !EQ(V[1], rhs.V[1], Precision) || !EQ(V[2], rhs.V[2], Precision));
}

/*!
*	Subscript operator for non-const object
*
*	\param idx index
*/
double &GVector3::operator [](const int &idx)
{
	assert(idx >= 0 && idx < 3);
	return V[idx];
}

/*!
*	Subscript operator for const object
*
*	\param idx index
*/
const double &GVector3::operator [](const int &idx) const
{
	assert(idx >= 0 && idx < 3);
	return V[idx];
}

/*!
*	Set the coordinates of a vector
*
*	\param x x-coordinate
*	\param y y-coordinate
*	\param z z-coordinate
*/
GVector3 &GVector3::Set(const double &x, const double &y, const double &z)
{
	V[0] = x;
	V[1] = y;
	V[2] = z;
	return *this;
}

/*!
*	Normalize the vector
*
*	\return reference to the vector
*/
GVector3 &GVector3::Normalize()
{
	double len = norm(*this);
	if (EQ_ZERO(len, Precision))
		return *this;
	V[0] /= len;
	V[1] /= len;
	V[2] /= len;
	return *this;
}

/*!
*	Set precision for equality and inequality test
*
*	\param error machine epsilon
*/
void GVector3::SetPrecision(double error)
{
	Precision = error;
}

/*!
*	Get current precision
*
*	\return current precision
*/
double GVector3::GetPrecision()
{
	return Precision;
}

/*!
*	Operator * (scaling)overloading
*
*	\param lhs left hand operand
*	\param s right hand operand
*	\return the vector scaled by \a s
*/
GVector3 operator *(const GVector3 &lhs, const double &s)
{
	GVector3 ret(lhs);
	ret *= s;
	return ret;
}

/*!
*	Operator * (scaling)overloading
*
*	\param s left hand operand
*	\param rhs right hand operand
*	\return the vector scaled by \a s
*/
GVector3 operator *(const double &s, const GVector3 &rhs)
{
	GVector3 ret(rhs);
	ret *= s;
	return ret;
}

/*!
*	Output operator
*
*	\param os output stream
*	\param v the output vector
*/
ostream &operator <<(ostream &os, const GVector3 &v)
{
	os << "(" << setw(10) << v.V[0] << ", " << setw(10) << v.V[1] << ", " << setw(10) << v.V[2] << ")";
	return os;
}

/*!
*	Input operator
*
*	\param is input stream
*	\param v the input vector
*/
istream &operator >> (istream &is, GVector3 &v)
{
	is >> v.V[0] >> v.V[1] >> v.V[2];
	return is;
}

/*!
*	\brief Project v to w.
*
*	\param v vector to project.
*	\param w target vector.
*
*	\return projection of v on w.
*/
GVector3 proj(const GVector3 &v, const GVector3 &w)
{
	return (v * w / SQR(norm(w))) * w;
}

/*!
*	\brief the distance from v to w.
*
*	\param v vector.
*	\param w vector.
*
*	\return the distance from v to w.
*/
double dist(const GVector3 &v, const GVector3 &w)
{
	return SQRT(SQR(v.V[0] - w.V[0]) + SQR(v.V[1] - w.V[1]) + SQR(v.V[2] - w.V[2]));
}

/*!
*	\brief size of the vector v.
*
*	\param v vector.
*
*	\return size of the vector v.
*/
double norm(const GVector3 &v)
{
	return SQRT(SQR(v.V[0]) + SQR(v.V[1]) + SQR(v.V[2]));
}

/*!
*	\brief angle between vector v and w.
*
*	\param v vector.
*	\param w vector.
*	\param radian true: radian, false: degree.
*
*	\return angle between vector v and w.
*/

double angle(const GVector3 &v, const GVector3 &w, bool radian)
{
	GVector3 p(v);
	GVector3 q(w);
	double cs, sn, angle;

	p.Normalize();
	q.Normalize();

	cs = p * q;
	sn = norm(p ^ q);

	angle = radian ? atan2(sn, cs) : RAD2DEG(atan2(sn, cs));
	return angle;
}