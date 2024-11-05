// Copyright Zhangci 2018
#include "GNURBS.h"
#include "NurbsProject.h"


double GKnot::Precision = 0.00000001;
double GNurbsCrv::Precision = 0.00000001;

/*!
*	\brief	default constructor
*
*	\param _m the last index of knot
*	\param _knot, knot = {x0, x1, ... , x_m}
*
*	\note default knot = {0, 0, 0, 0, 1, 1, 1, 1}.
*/
 
GKnot::GKnot(const int _m, const double *_knot)
{
	m = _m;
	knot = new double[m + 1];
	if (_knot != NULL)
		ARR_COPY(knot, _knot, m + 1);
	else
		ARR_ZERO(knot, m + 1);
}

 
GKnot::GKnot(const GKnot &copy)
{
	m = copy.m;
	knot = new double[m + 1];
	ARR_COPY(knot, copy.knot, m + 1);
}

 
GKnot::~GKnot()
{
	if (knot != NULL)
		delete[] knot;
	knot = NULL;
}

 
GKnot &GKnot::operator =(const GKnot &rhs)
{
	m = rhs.m;
	if (knot != NULL)
		delete[] knot;
	knot = new double[m + 1];
	ARR_COPY(knot, rhs.knot, m + 1);
	return *this;
}

 
double &GKnot::operator [](const int &idx)
{
	assert(idx >= 0 && idx <= m);
	return knot[idx];
}

/*!
*	\brief get the last inded of the knot.
*
*/
 
int GKnot::GetKnotIdx() const
{
	return m;
}

/*!
*	\brief get the knot vectors.
*/
 
double *GKnot::GetKnotVector()
{
	return knot;
}

/*!
*	\brief Different knot value of knot[st] between knot[ed].
*	\note Include knot[st] and knot[ed].
*
*	\param st 1st index.
*	\param ed the last index.
*	\param different idx.
*
*	\return knot vector.
*/
 
double *GKnot::GetSubVector(const int st, const int ed, int &idx) const
{
	static double subknot[500];
	int i;
	idx = 0;

	subknot[0] = knot[st];
	for (i = st + 1; i < ed + 1; i++)
	{
		if (knot[i] != subknot[idx])
			subknot[++idx] = knot[i];
	}

	return subknot;
}

 
const double &GKnot::operator [](const int &idx) const
{
	assert(idx >= 0 && idx <= m);
	return knot[idx];
}

/*!
*	\brief	knotspan
*/
 
int GKnot::FindKnotSpan(const int p, const double u) const
{
	assert(u >= knot[p] && u <= knot[m - p]);

	int low, high, mid, n;
	n = m - p - 1;
	if (u == knot[n + 1])
		return n;

	low = p;
	high = n + 1;
	mid = (low + high) / 2;
	while (u < knot[mid] || u >= knot[mid + 1])
	{
		if (u < knot[mid])
			high = mid;
		else
			low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

/*!
*	\brief	get multiplicity of u from Knot.
*
*	\param	dimension of the nurbs.
*	\param	u parameter.
*
*	\return multiplicity of u.
*/
 
int GKnot::FindKnotMult(const int p, const double u) const
{
	int mult = 0;
	for (int i = 0; i < m + 1; i++)
	{
		if (EQ(u, knot[i], Precision))
			mult++;
	}
	return mult;
}

/*!
*	\brief	insert r times by using u into the current knot.
*
*	\param p dimension of the nurbs using the knot.
*	\param u new vaule.
*	\param r multiplicity of insertion.
*	\param idx to save start index.
*	\param mul to save multiplicity before insertion.
*/
 
void GKnot::InsertKnots(const int p, const double u, int &r, int &idx, int &mul)
{
	int i;
	idx = FindKnotSpan(p, u);
	mul = FindKnotMult(p, u);

	if (u == knot[m - p])
		idx++;

	if (r + mul > p)
		r = p - mul;

	double *newKnot = new double[m + r + 1];

	// load new knot vectors.
	for (i = 0; i <= idx; i++)
		newKnot[i] = knot[i];

	for (i = idx + 1; i <= idx + r; i++)
		newKnot[i] = u;

	for (i = idx + r + 1; i <= m + r; i++)
		newKnot[i] = knot[i - r];

	delete[] knot;
	knot = newKnot;
	m = m + r;
}

 
int GKnot::RemoveKnots(const int p, const double u, int &r, GVector3 *P)
{
	assert(r >= 0);
	int idx = FindKnotSpan(p, u);
	int mul = FindKnotMult(p, u);


	if (mul <= 0)
		return 0;
	if (idx <= p || idx >= m - p)
		return 0;

	r = MIN(r, mul);

	int n = m - p - 1;
	int ord = p + 1;
	int fout = (2 * idx - mul - p) / 2;
	int last = idx - mul;
	int first = idx - p;

	int t, i, j, ii, jj, off, remflag;
	double alfi, alfj;
	int k;

	GVector3 *tmp = new GVector3[2 * p + 1];

	for (t = 0; t < r; t++)
	{
		off = first - 1;
		tmp[0] = P[off];
		tmp[last + 1 - off] = P[last + 1];

		i = first;
		j = last;
		ii = 1;
		jj = last - off;
		remflag = 0;

		while (j - i > t)
		{
			alfi = (u - knot[i]) / (knot[i + ord + t] - knot[i]);
			alfj = (u - knot[j - t]) / (knot[j + ord] - knot[j - t]);

			tmp[ii] = (P[i] - (1.0 - alfi) * tmp[ii - 1]) / alfi;
			tmp[jj] = (P[j] - alfj * tmp[jj + 1]) / (1.0 - alfj);

			i++; ii++;
			j--; jj--;
		}

		if (j - i < t)
		{
			if (tmp[ii - 1] == tmp[jj + 1])
				remflag = 1;
		}
		else
		{
			alfi = (u - knot[i]) / (knot[i + ord + t] - knot[i]);
			GVector3 X;
			X = tmp[ii - 1] * (1.0 - alfi) + tmp[ii + t + 1] * alfi;
			if (P[i] == X)
				remflag = 1;
		}

		if (remflag == 0)
			break;
		else
		{
			i = first;
			j = last;
			while (j - i > t)
			{
				P[i] = tmp[i - off];
				P[j] = tmp[j - off];
				i++; j--;
			}
		}
		first--;
		last++;
	}

	if (t == 0)
	{
		delete[] tmp;
		return t;
	}

	m = m - t;
	double *newKnot = new double[m + 1];

	for (k = 0; k < idx - t + 1; k++)
		newKnot[k] = knot[k];
	for (k = idx + 1; k <= m + t; k++)
		newKnot[k - t] = knot[k];

	delete[] knot;
	knot = newKnot;

	j = fout;
	i = j;
	for (k = 1; k < t; k++)
	{
		if (k % 2 == 1)
			i++;
		else
			j--;
	}

	for (k = i + 1; k <= n; k++)
	{
		P[j] = P[k];
		j++;
	}
	delete[] tmp;
	return t;
}

/*!
*	\brief Compute basis of nurbs.
*
*	\param p dimension of the nurbs.
*	\param u
*	\param idx start index of frist parameter u.
*	\param basis
*/

void GKnot::GetBasis(const int p, const double u, const int idx, double *basis) const
{
	int j, k;
	basis[0] = 1.0;
	for (j = 1; j <= p; j++)
	{
		double saved = 0.0;
		for (k = 0; k < j; k++)
		{
			double tmp = basis[k] / (knot[idx + k + 1] - knot[idx + k + 1 - j]);
			basis[k] = saved + (knot[idx + k + 1] - u) * tmp;
			saved = (u - knot[idx + 1 - j + k]) * tmp;
		}
		basis[j] = saved;
	}
}

/*!
*	\brief	Compute Derivative value of basis

*/
 
void GKnot::GetDerivBasis(const int p, const double u, const int idx, const int nth, double *basis) const
{
	static double drvs[MAX_DEG][MAX_DEG];
	static double ndu[MAX_DEG][MAX_DEG];
	static double left[MAX_DEG];
	static double right[MAX_DEG];
	static double a[2][MAX_DEG];

	double saved, tmp, dd;
	int j, r, k, s1, s2, j1, j2, rk, pk;
	ndu[0][0] = 1.0;
	for (j = 1; j < p + 1; j++)
	{
		left[j] = u - knot[idx + 1 - j];
		right[j] = knot[idx + j] - u;
		saved = 0.0;
		for (r = 0; r < j; r++)
		{
			ndu[j][r] = right[r + 1] + left[j - r];
			tmp = ndu[r][j - 1] / ndu[j][r];
			ndu[r][j] = saved + right[r + 1] * tmp;
			saved = left[j - r] * tmp;
		}
		ndu[j][j] = saved;
	}

	for (j = 0; j < p + 1; j++)
		drvs[0][j] = ndu[j][p];

	for (r = 0; r < p + 1; r++)
	{
		s1 = 0;
		s2 = 1;
		a[0][0] = 1.0;
		for (k = 1; k < p + 1; k++)
		{
			dd = 0.0;
			rk = r - k, pk = p - k;
			if (r >= k)
			{
				a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
				dd = a[s2][0] * ndu[rk][pk];
			}
			if (rk >= -1)
				j1 = 1;
			else
				j1 = -rk;
			if (r - 1 <= pk)
				j2 = k - 1;
			else
				j2 = p - r;
			for (j = j1; j <= j2; j++)
			{
				a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
				dd += a[s2][j] * ndu[rk + j][pk];
			}
			if (r <= pk)
			{
				a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
				dd += a[s2][k] * ndu[r][pk];
			}
			drvs[k][r] = dd;
			j = s1; s1 = s2; s2 = j;
		}
	}
	r = p;
	for (k = 1; k < p + 1; k++)
	{
		for (j = 0; j < p + 1; j++)
		{
			drvs[k][j] *= r;
		}
		r *= (p - k);
	}

	ARR_COPY(basis, drvs[nth], p + 1);
}


 
void GKnot::SetPrecision(double error)
{
	Precision = error;
}


 
double GKnot::GetPrecision()
{
	return Precision;
}



/*!
*/
 
GNurbsCrv::GNurbsCrv()
{
	p = 0;
	n = 0;
	P = new GVector3[1];
	P[0].Set(0.0, 0.0, 0.0);

	double knot[2] = { 0.0, 1.0 };
	U = GKnot(n + p + 1, knot);
	bClosed = false;
}

/*!
*	\brief	create nurbs curve.
*
*	\param p degree
*	\param n the last index of control points.
*	\param P control points
*	\param U knot vector
*
*	\note the dimension of control point can be arbitrary to create different curves.
*/
 
GNurbsCrv::GNurbsCrv(const int _p, const int _n, const GVector3 *_P, const double *_knot)
	: U(_n + _p + 1, _knot)
{
	p = _p;
	n = _n;
	P = new GVector3[n + 1];
	for (int i = 0; i < n + 1; i++)
		P[i] = _P[i];

	bClosed = false;
}


 
GNurbsCrv::GNurbsCrv(const GNurbsCrv &copy)
	: U(copy.U)
{
	p = copy.p;
	n = copy.n;

	P = new GVector3[n + 1];
	for (int i = 0; i < n + 1; i++)
		P[i] = copy.P[i];

	bClosed = copy.bClosed;
}

/*!
*
*/
 
GNurbsCrv::~GNurbsCrv()
{
	if (P != NULL)
		delete[] P;
	P = NULL;
}


 
GNurbsCrv &GNurbsCrv::operator =(const GNurbsCrv &rhs)
{
	p = rhs.p;
	n = rhs.n;

	if (P != NULL)
		delete[] P;
	P = new GVector3[n + 1];
	for (int i = 0; i < n + 1; i++)
		P[i] = rhs.P[i];

	U = rhs.U;
	bClosed = rhs.bClosed;

	return *this;
}

/*!
*	\brief	compute point on nurbs curve.
*
*	\param u
*	\param nth derivative order
*
*	\return vector
*/
 
GVector3 GNurbsCrv::Eval(const double &u, const int &nth) const
{
	int idx = U.FindKnotSpan(p, u);
	static double basis[MAX_DEG];

	if (nth == 0)
		U.GetBasis(p, u, idx, basis);
	else
		U.GetDerivBasis(p, u, idx, nth, basis);

	int dim = GetDim();
	GVector3 C(dim);
	for (int i = 0; i < p + 1; i++)
		C += basis[i] * P[idx - p + i];

	return C;
}

 
int GNurbsCrv::GetDegree() const
{
	return p;
}

 
int GNurbsCrv::GetDim() const
{
	//return P[0].GetDim();
	return 3;
}


 
int GNurbsCrv::GetCtlPtIdx() const
{
	return n;
}


 
GVector3 *GNurbsCrv::GetCtlPt() const
{
	return P;
}

 
GKnot GNurbsCrv::GetKnots() const
{
	return U;
}

 
void GNurbsCrv::GetDomain(double &min, double &max) const
{
	min = U[p];
	max = U[n + 1];
}

 
bool GNurbsCrv::IsClosed() const
{
	return bClosed;
}


 
bool GNurbsCrv::IsBezierForm() const
{
	return (p == n);
}


 
void GNurbsCrv::InsertKnots(const double u, int r)
{
	assert(r > 0);
	int i, j, L = 0;
	int idx;			// knot span index.
	int s;				// knot multiplicity.
	GKnot oldU = U;

	// Update knot vector.
	// r is decreased if the sum of r and u's multiplicity is greater than p.
	U.InsertKnots(p, u, r, idx, s);
	if (r <= 0)
		return;

	GVector3 *P1 = new GVector3[n + r + 1];		// new control points.
	GVector3 *R = new GVector3[p + 1];			// temporary control points.

												// save unaltered control points
	for (i = 0; i <= idx - p; i++)
		P1[i] = P[i];

	for (i = idx - s; i <= n; i++)
		P1[i + r] = P[i];

	for (i = 0; i <= p - s; i++)
		R[i] = P[idx - p + i];

	// knot insert process.
	for (j = 1; j <= r; j++)
	{
		L = idx - p + j;
		for (i = 0; i <= p - j - s; i++)
		{
			double alpha = (u - oldU[L + i]) / (oldU[i + idx + 1] - oldU[L + i]);
			R[i] = alpha * R[i + 1] + (1.0 - alpha) * R[i];
		}
		P1[L] = R[0];
		P1[idx + r - j - s] = R[p - j - s];
	}
	for (i = L + 1; i < idx - s; i++)
		P1[i] = R[i - L];

	delete[] P;
	delete[] R;
	P = P1;
	n = n + r;
}

/*!
*	\brief	Remove knot from knot vector.
*
*	\param u kont to remove.
*	\param r multiplicity of removal.
*
*	\return multiplicity of removal.
*/
 
int GNurbsCrv::RemoveKnots(const double u, int r)
{
	int scs = U.RemoveKnots(p, u, r, P);

	if (scs == 0)
		return scs;

	int new_n = n - scs;

	GVector3 *nP = new GVector3[new_n + 1];
	for (int i = 0; i < new_n + 1; i++)
		nP[i] = P[i];
	delete[] P;

	n = new_n;
	P = nP;

	return scs;
}

/*!
*	\brief refine knot vector
*/
 
void GNurbsCrv::Refinement()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p, n + 1, idx);

	for (i = 0; i < idx; i++)
	{
		double u = (knot[i] + knot[i + 1]) * 0.5;
		InsertKnots(u, 1);
	}
}

 
void GNurbsCrv::MakeBezierForm()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p + 1, n, idx);

	for (i = 0; i < idx + 1; i++)
		InsertKnots(knot[i], p);
}


 
void GNurbsCrv::MakeCompactForm()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p + 1, n, idx);
	for (i = 0; i < idx + 1; i++)
		RemoveKnots(knot[i], p);
}

 
void GNurbsCrv::SetPrecision(double error)
{
	Precision = error;
}

 
double GNurbsCrv::GetPrecision()
{
	return Precision;
}

/*!
*	\brief Edit a point of curve by gaving an offset.
*
*	\param u curve parameter C(u)=...
*	\param offset
*/
 
void GNurbsCrv::Edit(const double &u, const GVector3 &offset)
{
	int i, j;

	// get nodes.
	double *node = new double[n + 1];
	ARR_ZERO(node, n + 1);

	for (i = 0; i < n + 1; i++)
	{
		double sum = 0.0;
		for (j = 1; j < p + 1; j++)
			sum += U[i + j];
		node[i] = sum / p;
	}

	// get index k, P[k] is translated to produce the desired movement.
	int k = 0;
	double min = 10000000.0f;
	for (i = 0; i < n + 1; i++)
	{
		double tmp = fabs(u - node[i]);
		if (tmp < min)
		{
			k = i;
			min = tmp;
		}
	}

	int idx = U.FindKnotSpan(p, u);
	static double basis[MAX_DEG];
	U.GetBasis(p, u, idx, basis);

	double Rmax = -10000.0;
	for (i = 0; i < p + 1; i++)
		Rmax = MAX(basis[i], Rmax);

	Rmax = 1 / Rmax;
	P[k] = P[k] + Rmax * offset;

	if (bClosed)
	{
		if (k >= 0 && k <= p - 1)
			P[k + n - p + 1] = P[k + n - p + 1] + Rmax * offset;
		if (k >= n - p + 1 && k <= n)
			P[k - n + p - 1] = P[k - n + p - 1] + Rmax * offset;
	}

	// memory free
	delete[] node;
}


GNurbsCrv *get_gnurbs_closed_crv(int _p, int _n, GVector3 *_P)
{
	int p = _p;
	int n = _n + p;
	int m = n + p + 1;

	GVector3 *P = new GVector3[n + 1];
	for (int i = 0; i < n + 1; ++i)
	{
		int idx = i % (_n + 1);
		P[i] = _P[idx];
	}

	double *U = new double[m + 1];
	for (int i = 0; i < m + 1; ++i)
		U[i] = (double)i;

	double a = U[p], b = U[n + 1];
	for (int i = 0; i < m + 1; ++i)
		U[i] = (U[i] - a) / (b - a);

	GNurbsCrv *pCrv = new GNurbsCrv(p, n, P, U);
	pCrv->bClosed = true;

	delete[] P;
	delete[] U;

	return pCrv;
}

GNurbsSrf::GNurbsSrf()
{
	double knot[2] = { 0.0, 1.0 };

	p = 0;
	m = 0;
	U = GKnot(m + p + 1, knot);

	q = 0;
	n = 0;
	V = GKnot(n + q + 1, knot);

	P = new GVector3 *[1];
	P[0] = new GVector3[1];
	P[0][0].Set(0.0, 0.0, 0.0);
	bClosedU = bClosedV = false;
}

GNurbsSrf::GNurbsSrf(const int _p, const int _m, const double *_U, const int _q, const int _n, const double *_V, const GVector3 *_P)
	: U(_p + _m + 1, _U), V(_q + _n + 1, _V)
{
	p = _p;
	m = _m;
	q = _q;
	n = _n;
	P = new GVector3 *[m + 1];
	for (int i = 0; i < m + 1; i++)
	{
		P[i] = new GVector3[n + 1];
		for (int j = 0; j < n + 1; j++)
			P[i][j] = _P[(n + 1) * i + j];
	}
	bClosedU = bClosedV = false;
}

GNurbsSrf::GNurbsSrf(const GNurbsSrf &copy)
	: U(copy.U), V(copy.V)
{
	p = copy.p;
	m = copy.m;
	q = copy.q;
	n = copy.n;
	P = new GVector3 *[m + 1];
	for (int i = 0; i < m + 1; i++)
	{
		P[i] = new GVector3[n + 1];
		for (int j = 0; j < n + 1; j++)
			P[i][j] = copy.P[i][j];
	}
	bClosedU = copy.bClosedU;
	bClosedV = copy.bClosedV;
}


GNurbsSrf::~GNurbsSrf()
{
	if (P != NULL)
	{
		for (int i = 0; i < m + 1; i++)
			delete[] P[i];
		delete[] P;
		P = NULL;
	}
}

GNurbsSrf &GNurbsSrf::operator =(const GNurbsSrf &rhs)
{
	p = rhs.p;
	m = rhs.m;
	q = rhs.q;
	n = rhs.n;
	U = rhs.U;
	V = rhs.V;
	if (P != NULL)
	{
		for (int i = 0; i < m + 1; i++)
			delete[] P[i];
		delete[] P;
	}
	P = new GVector3 *[m + 1];
	for (int i = 0; i < m + 1; i++)
	{
		P[i] = new GVector3[n + 1];
		for (int j = 0; j < n + 1; j++)
			P[i][j] = rhs.P[i][j];
	}

	bClosedU = rhs.bClosedU;
	bClosedV = rhs.bClosedV;

	return *this;
}

GVector3 GNurbsSrf::Eval(const double &u, const double &v, const int &kth, const int &lth) const
{
	static double uBasis[MAX_DEG];
	static double vBasis[MAX_DEG];
	int uIdx = U.FindKnotSpan(p, u);
	int vIdx = V.FindKnotSpan(q, v);

	if (kth == 0)
		U.GetBasis(p, u, uIdx, uBasis);
	else
		U.GetDerivBasis(p, u, uIdx, kth, uBasis);

	if (lth == 0)
		V.GetBasis(q, v, vIdx, vBasis);
	else
		V.GetDerivBasis(q, v, vIdx, lth, vBasis);

	int dim = GetDim();

	GVector3 S(dim);

	
	for (int i = 0; i < p + 1; i++)
	{
		GVector3 tmp1(dim);
		for (int j = 0; j < q + 1; j++)
		{
			int s = uIdx - p + i;
			int t = vIdx - q + j;
			tmp1 += vBasis[j] * P[s][t];
		}
		S += uBasis[i] * tmp1;
	}
	
	return S;
}

void GNurbsSrf::InsertKnotsU(const double u, int r)
{
	assert(r > 0);
	int i, j, k, L=0;		// index 
	int idx;			// knot span index.
	int mul;			// knot multiplicity.
	GKnot oldU = U;		// save old knot vector.

						// Update knot vector.
	U.InsertKnots(p, u, r, idx, mul);

	GVector3 **P1 = new GVector3 *[m + r + 1];		// new control points.
	for (i = 0; i < m + r + 1; i++)
		P1[i] = new GVector3[n + 1];

	GVector3 *R = new GVector3[p + 1];			// temporary control points.

	double **alpha = new double *[p];
	for (i = 0; i < p; i++)
		alpha[i] = new double[r + 1];

	// compute the alphas
	for (j = 1; j <= r; j++)
	{
		L = idx - p + j;
		for (i = 0; i <= p - j - mul; i++)
			alpha[i][j] = (u - oldU[L + i]) / (oldU[i + idx + 1] - oldU[L + i]);
	}

	for (k = 0; k < n + 1; k++)
	{
		// save unaltered control points
		for (i = 0; i <= idx - p; i++)
			P1[i][k] = P[i][k];

		for (i = idx - mul; i <= m; i++)
			P1[i + r][k] = P[i][k];

		for (i = 0; i <= p - mul; i++)
			R[i] = P[idx - p + i][k];

		// knot insert process.
		for (j = 1; j <= r; j++)
		{
			L = idx - p + j;
			for (i = 0; i <= p - j - mul; i++)
				R[i] = alpha[i][j] * R[i + 1] + (1.0 - alpha[i][j]) * R[i];
			P1[L][k] = R[0];
			P1[idx + r - j - mul][k] = R[p - j - mul];
		}

		for (i = L + 1; i < idx - mul; i++)
			P1[i][k] = R[i - L];
	}

	// delete alpha
	for (i = 0; i < p; i++)
		delete[] alpha[i];
	delete[] alpha;

	// delete R
	delete[] R;

	// delete old control points P
	for (i = 0; i < m + 1; i++)
		delete[] P[i];
	delete[] P;

	// Update the control points and its number.
	P = P1;
	m = m + r;
}

void GNurbsSrf::InsertKnotsV(const double v, int r)
{
	assert(r > 0);
	int i, j, k, L=0;		// index 
	int idx;			// knot span index.
	int mul;			// knot multiplicity.
	GKnot oldV = V;		// save old knot vector.

						// Update knot vector.
	V.InsertKnots(q, v, r, idx, mul);

	GVector3 **P1 = new GVector3 *[m + 1];		// new control points.
	for (i = 0; i < m + 1; i++)
		P1[i] = new GVector3[n + r + 1];

	GVector3 *R = new GVector3[q + 1];			// temporary control points.

	double **alpha = new double *[q];
	for (i = 0; i < q; i++)
		alpha[i] = new double[r + 1];

	// compute the alphas
	for (j = 1; j <= r; j++)
	{
		L = idx - q + j;
		for (i = 0; i <= q - j - mul; i++)
			alpha[i][j] = (v - oldV[L + i]) / (oldV[i + idx + 1] - oldV[L + i]);
	}

	for (k = 0; k < m + 1; k++)
	{
		// save unaltered control points
		for (i = 0; i <= idx - q; i++)
			P1[k][i] = P[k][i];

		for (i = idx - mul; i <= n; i++)
			P1[k][i + r] = P[k][i];

		for (i = 0; i <= q - mul; i++)
			R[i] = P[k][idx - q + i];

		// knot insert process.
		for (j = 1; j <= r; j++)
		{
			L = idx - q + j;
			for (i = 0; i <= q - j - mul; i++)
				R[i] = alpha[i][j] * R[i + 1] + (1.0 - alpha[i][j]) * R[i];
			P1[k][L] = R[0];
			P1[k][idx + r - j - mul] = R[q - j - mul];
		}

		for (i = L + 1; i < idx - mul; i++)
			P1[k][i] = R[i - L];
	}

	// delete alpha
	for (i = 0; i < p; i++)
		delete[] alpha[i];
	delete[] alpha;

	// delete R
	delete[] R;

	// delete old control points P
	for (i = 0; i < m + 1; i++)
		delete[] P[i];
	delete[] P;

	// Update the control points and its number.
	P = P1;
	n = n + r;
}

int GNurbsSrf::RemoveKnotsU(const double u, int r)
{
	assert(r > 0);

	int i, j, k, tot_scs = 0;
	GVector3 **Q = new GVector3 *[n + 1];
	for (i = 0; i < n + 1; i++)
	{
		Q[i] = new GVector3[m + 1];
		for (j = 0; j < m + 1; j++)
			Q[i][j] = P[j][i];
	}

	for (i = 0; i < r; i++)
	{
		int scs = 0;
		GKnot tmpU;

		for (j = 0; j < n + 1; j++)
		{
			tmpU = U;

			int rnum = 1;
			scs = tmpU.RemoveKnots(p, u, rnum, Q[j]);

			if (scs != 1)
				break;
		}

		if (scs != 1)
			break;

		tot_scs++;

		for (j = 0; j < m + 1; j++)
			delete[] P[j];
		delete[] P;

		m = m - 1;

		P = new GVector3 *[m + 1];
		for (j = 0; j < m + 1; j++)
			P[j] = new GVector3[n + 1];

		for (j = 0; j < m + 1; j++)
			for (k = 0; k < n + 1; k++)
				P[j][k] = Q[k][j];

		U = tmpU;
	}

	for (i = 0; i < n + 1; i++)
		delete[] Q[i];
	delete[] Q;

	return tot_scs;
}

int GNurbsSrf::RemoveKnotsV(const double v, int r)
{
	assert(r > 0);

	int i, j, k, tot_scs = 0;
	GVector3 **Q = new GVector3 *[m + 1];
	for (i = 0; i < m + 1; i++)
	{
		Q[i] = new GVector3[n + 1];
		for (j = 0; j < n + 1; j++)
			Q[i][j] = P[i][j];
	}

	for (i = 0; i < r; i++)
	{
		int scs = 0;
		GKnot tmpV;

		for (j = 0; j < m + 1; j++)
		{
			tmpV = V;

			int rnum = 1;
			scs = tmpV.RemoveKnots(q, v, rnum, Q[j]);

			if (scs != 1)
				break;
		}

		if (scs != 1)
			break;

		tot_scs++;

		for (j = 0; j < m + 1; j++)
			delete[] P[j];
		delete[] P;

		n = n - 1;

		P = new GVector3 *[m + 1];
		for (j = 0; j < m + 1; j++)
			P[j] = new GVector3[n + 1];

		for (j = 0; j < m + 1; j++)
			for (k = 0; k < n + 1; k++)
				P[j][k] = Q[j][k];

		V = tmpV;
	}

	for (i = 0; i < m + 1; i++)
		delete[] Q[i];
	delete[] Q;

	return tot_scs;
}

void GNurbsSrf::RefinementU()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p, m + 1, idx);

	for (i = 0; i < idx; i++)
	{
		double u = (knot[i] + knot[i + 1]) * 0.5;
		InsertKnotsU(u, 1);
	}
}

void GNurbsSrf::RefinementV()
{
	int i, idx = 0;
	double *knot;

	knot = V.GetSubVector(q, n + 1, idx);

	for (i = 0; i < idx; i++)
	{
		double v = (knot[i] + knot[i + 1]) * 0.5;
		InsertKnotsV(v, 1);
	}
}

void GNurbsSrf::MakeBezierForm()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p + 1, m, idx);
	for (i = 0; i < idx + 1; i++)
		InsertKnotsU(knot[i], p);

	knot = V.GetSubVector(q + 1, n, idx);
	for (i = 0; i < idx + 1; i++)
		InsertKnotsV(knot[i], q);
}

void GNurbsSrf::MakeCompactForm()
{
	int i, idx = 0;
	double *knot;

	knot = U.GetSubVector(p + 1, m, idx);
	for (i = 0; i < idx + 1; i++)
		RemoveKnotsU(knot[i], p);

	knot = V.GetSubVector(q + 1, n, idx);
	for (i = 0; i < idx + 1; i++)
		RemoveKnotsV(knot[i], q);
}

void GNurbsSrf::Edit(const double &u, const double &v, const GVector3 &offset)
{
	int i, j;

	// get nodes.
	double *node_u = new double[m + 1];
	double *node_v = new double[n + 1];
	ARR_ZERO(node_u, m + 1);
	ARR_ZERO(node_v, n + 1);

	for (i = 0; i < m + 1; i++)
	{
		double sum = 0.0;
		for (j = 1; j < p + 1; j++)
			sum += U[i + j];
		node_u[i] = sum / p;
	}

	for (i = 0; i < n + 1; i++)
	{
		double sum = 0.0;
		for (j = 1; j < q + 1; j++)
			sum += V[i + j];
		node_v[i] = sum / q;
	}

	// get index k,l, P[k][l] is translated to produce the desired movement.
	int k = 0, l = 0;
	double min = 10000000.0f;
	for (i = 0; i < m + 1; i++)
	{
		double tmp = fabs(u - node_u[i]);
		if (tmp < min)
		{
			k = i;
			min = tmp;
		}
	}

	min = 1000000.0f;
	for (i = 0; i < n + 1; i++)
	{
		double tmp = fabs(v - node_v[i]);
		if (tmp < min)
		{
			l = i;
			min = tmp;
		}
	}

	int uIdx = U.FindKnotSpan(p, u);
	int vIdx = V.FindKnotSpan(q, v);
	static double uBasis[MAX_DEG], vBasis[MAX_DEG];

	U.GetBasis(p, u, uIdx, uBasis);
	V.GetBasis(q, v, vIdx, vBasis);

	double Rmax_U = -10000.0;
	double Rmax_V = -10000.0;
	double w = 0.0;

	for (i = 0; i < p + 1; i++)
		Rmax_U = MAX(uBasis[i], Rmax_U);

	for (i = 0; i < q + 1; i++)
		Rmax_V = MAX(vBasis[i], Rmax_V);


	for (i = 0; i < q + 1; i++)
	{
		double tmp = 0.0;
		for (j = 0; j < p + 1; j++)
		{
			tmp += uBasis[j];
		}
		w += vBasis[i] * tmp; //w
	}
	double alpha = Rmax_U * Rmax_V / w;
	alpha = 1.0 / alpha;

	P[k][l] += alpha * offset;

	if (bClosedU)
	{
		if (k >= 0 && k <= p - 1)
			P[k + m - p + 1][l] += alpha * offset;
		if (k >= m - p + 1 && k <= m)
			P[k - m + p - 1][l] += alpha * offset;
	}

	if (bClosedV)
	{
		if (l >= 0 && l <= q - 1)
			P[k][l + n - q + 1] += alpha * offset;
		if (l >= n - q + 1 && l <= n)
			P[k][l - n + q - 1] += alpha * offset;
	}

	// memory free
	delete[] node_u;
	delete[] node_v;
}

GNurbsSrf *get_gnurbs_closed_srf(int _p, int _m, int _q, int _n, GVector3 *_P, bool _bClosedU, bool _bClosedV)
{
	int p = _p;
	int m = (_bClosedU) ? _m + p : _m;
	int r = m + p + 1;
	int q = _q;
	int n = (_bClosedV) ? _n + q : _n;
	int s = n + q + 1;

	GVector3 *P = new GVector3[m * n + m + n + 1];

	for (int i = 0; i < m + 1; ++i)
	{
		int idx0, idx1;
		idx0 = (_bClosedU) ? i % (_m + 1) : i;

		for (int j = 0; j < n + 1; ++j)
		{
			idx1 = (_bClosedV) ? j % (_n + 1) : j;
			P[i * (n + 1) + j] = _P[idx0 * (_n + 1) + idx1];
		}
	}

	double a, b;
	double *U = new double[r + 1];
	double *V = new double[s + 1];

	for (int i = 0; i < r + 1; ++i)
		U[i] = (double)i;

	for (int i = 0; i < s + 1; ++i)
		V[i] = (double)i;

	a = U[p];
	b = U[m + 1];
	for (int i = 0; i < r + 1; ++i)
		U[i] = (U[i] - a) / (b - a);

	a = V[q];
	b = V[n + 1];
	for (int i = 0; i < s + 1; ++i)
		V[i] = (V[i] - a) / (b - a);

	GNurbsSrf *pSrf = new GNurbsSrf(p, m, U, q, n, V, P);
	pSrf->bClosedU = _bClosedU;
	pSrf->bClosedV = _bClosedV;

	delete[] P;
	delete[] U;
	delete[] V;

	return pSrf;
}