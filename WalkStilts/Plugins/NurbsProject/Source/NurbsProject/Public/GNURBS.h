// Copyright Zhangci 2018
//class of NURBS's algorithm
//brief Basic algorithm of NURBS, included knot, curve and surface.

#pragma once
#include <iostream>
#include "gmath.h"

#define MAX_DEG 10

class GKnot;
class GNurbsCrv;
class GNurbsSrf;

/*!
*	\class GKnot
*	\brief knot vector of nurbs

*/

class GKnot
{
public:
	GKnot(const int _m = 7, const double *_knot = NULL);
	GKnot(const GKnot &copy);
	virtual ~GKnot();

	GKnot &operator =(const GKnot &rhs);
	double &operator [](const int &idx);
	const double &operator [](const int &idx) const;

	int		GetKnotIdx() const;
	double *GetKnotVector();
	double *GetSubVector(const int st, const int ed, int &idx) const;

	int		FindKnotSpan(const int p, const double u) const;
	int		FindKnotMult(const int p, const double u) const;
	void	InsertKnots(const int p, const double u, int &r, int &idx, int &mul);
	int		RemoveKnots(const int p, const double u, int &r, GVector3 *P);
	void	GetBasis(const int p, const double u, const int idx, double *basis) const;
	void	GetDerivBasis(const int p, const double u, const int idx, const int nth, double *basis) const;

	static void SetPrecision(double error);
	static double GetPrecision();

	friend class GNurbsCrv;
	friend class GNurbsSrf;

protected:
	int m;				// count of knots: m + 1
	double *knot;		// knot vector, knot = {x0, x1, x2, ... , x_m}
	
	static double Precision;
};


/*!
*	\class GNurbsCrv
*	\brief Class of 3D Nurbs curves.
*/
class GNurbsCrv 
{
public:
	GNurbsCrv();
	GNurbsCrv(const int _p, const int _n, const GVector3 *_P, const double *_knot);
	GNurbsCrv(const GNurbsCrv &copy);
	virtual ~GNurbsCrv();

	GNurbsCrv &operator =(const GNurbsCrv &rhs);

	int		GetDim() const;
	int		GetDegree() const;
	int		GetCtlPtIdx() const;
	void	GetDomain(double &min, double &max) const;
	bool	IsClosed() const;
	bool	IsBezierForm() const;
	GKnot	GetKnots() const;
	GVector3 *GetCtlPt() const;

	void InsertKnots(const double u, int r);
	int	 RemoveKnots(const double u, int r);
	void Refinement();
	void MakeBezierForm();
	void MakeCompactForm();
	void Edit(const double &u, const GVector3 &offset);
	GVector3 Eval(const double &u, const int &nth = 0) const;

	static void SetPrecision(double error);
	static double GetPrecision();

	friend GNurbsCrv *get_gnurbs_closed_crv(int _p, int _n, GVector3 *_P);

protected:
	int p;				// degree of the nurbs curve
	int n;				// the last index of control points
	GVector3 *P;		// control points
	GKnot U;			// knot vector
	bool bClosed;		// Closed curve: true, Open curve: flase
	static double Precision;
};

GNurbsCrv *get_gnurbs_closed_crv(int _p, int _n, GVector3 *_P);


/*!
*	\class GNurbsSrf
*	\brief Class of 3D Nurbs surfaces.
*
*/
class GNurbsSrf
{
public:
	
	GNurbsSrf();

	/*!
	*	\brief	Create nurbs surface.
	*
	*	\param _p degree of u direction
	*	\param _m last index of control point for u direction
	*	\param _U knot of u direction
	*	\param _q degree of v direction
	*	\param _n last index of control point for v direction.
	*	\param _V knot of v direction
	*	\param _P control points, (m + 1) x (n + 1).
	*/
	GNurbsSrf(const int _p, const int _m, const double *_U, const int _q, const int _n, const double *_V, const GVector3 *_P);
	
	GNurbsSrf(const GNurbsSrf &copy);
	virtual ~GNurbsSrf();

	
	GNurbsSrf &operator =(const GNurbsSrf &rhs);


	int		GetDim() const;
	int		GetDegU() const;
	int		GetDegV() const;
	int		GetCtlPtIdxU() const;
	int		GetCtlPtIdxV() const;
	void	GetDomainU(double &min, double &max) const;
	void	GetDomainV(double &min, double &max) const;
	GKnot	GetKnotU() const;
	GKnot	GetKnotV() const;
	GVector3 **GetCtlPt() const;

	
	//void	SetQuatSrf();

	/*!
	*	\brief	Insert a knot into U direction knot vector.
	*
	*	\param u knot to insertion.
	*	\param r multiplicity of insertion.
	*/
	void	InsertKnotsU(const double u, int r = 1);

	/*!
	*	\brief	Insert a knot into V direction knot vector.
	*
	*	\param v knot to insertion.
	*	\param r multiplicity of insertion.
	*/
	void	InsertKnotsV(const double v, int r = 1);

	/*!
	*	\brief	Remove u direction knot from U direction knot vector.
	*
	*	\param u knot to removal.
	*	\param r multiplicity of removal.
	*
	*	\return actual multiplicity of removal.
	*/
	int		RemoveKnotsU(const double u, int r);

	/*!
	*	\brief	Remove v direction knot from V direction knot vector.
	*
	*	\param v knot to removal.
	*	\param r multiplicity of removal.
	*
	*	\return actual multiplicity of removal.
	*/
	int		RemoveKnotsV(const double v, int r);

	/*!
	*	\brief Insert middle value between all U direction knot vector.
	*/
	void	RefinementU();

	/*!
	*	\brief Insert middle value between all V direction knot vector.
	*/
	void	RefinementV();


	bool	IsBezierForm() const;
	bool	IsClosedU() const;
	bool	IsClosedV() const;
	void	MakeBezierForm();
	void	MakeCompactForm();

	/*!
	*	\brief find a proper control point to move surface point s(u,v).
	*
	*	\param u: u direction parameter.
	*	\param v: v direction parameter.
	*	\param offset: move offset.
	*/
	void	Edit(const double &u, const double &v, const GVector3 &offset);

	virtual GVector3 Eval(const double &u, const double &v, const int &kth = 0, const int &lth = 0) const;

	static void SetPrecision(double error);
	static double GetPrecision();


	/*!
	*	\brief Create close nurbs surface.
	*
	*	\param _p degree of u direction.
	*	\param _m control points last index of u direction. 
	*	\param _q degree of v direction. 
	*	\param _n control points last index of v direction.
	*	\param _P control points.
	*	\param _bClosedU u direction close
	*	\param _bClosedV v direction close.
	*
	*	\return close nurbs surface.
	*/
	friend GNurbsSrf *get_gnurbs_closed_srf(int _p, int _m, int _q, int _n, GVector3 *_P, bool _bClosedU, bool _bClosedV);

protected:
	
	int p;			// degree of u direction.
	int q;			// degree of v direction.
	int m;			// last index of control point for u direction. 
	int n;			// last index of control point for v direction. 
	GKnot U;		// knot of u direction.
	GKnot V;		// knot of v direction.
	GVector3 **P;	// control points, (m + 1) x (n + 1).
	bool bClosedU;	// closed surface of u direction: true, open surface: flase
	bool bClosedV;	// closed surface of v direction: true, open surface: flase
	//bool bQuatSrf;	// quaternion surface

	static double Precision;
};

inline
bool GNurbsSrf::IsBezierForm() const
{
	return ((p == m) && (q == n));
}


inline
bool GNurbsSrf::IsClosedU() const
{
	return bClosedU;
}


inline
bool GNurbsSrf::IsClosedV() const
{
	return bClosedV;
}

inline
int GNurbsSrf::GetDim() const
{
	return  3;
}

inline
int GNurbsSrf::GetDegU() const
{
	return p;
}

inline
int GNurbsSrf::GetDegV() const
{
	return q;
}

inline
GKnot GNurbsSrf::GetKnotU() const
{
	return U;
}

inline
GKnot GNurbsSrf::GetKnotV() const
{
	return V;
}

inline
int GNurbsSrf::GetCtlPtIdxU() const
{
	return m;
}

inline
int GNurbsSrf::GetCtlPtIdxV() const
{
	return n;
}

inline
GVector3 **GNurbsSrf::GetCtlPt() const
{
	return P;
}

inline
void GNurbsSrf::GetDomainU(double &min, double &max) const
{
	min = U[p];
	max = U[m + 1];
}

inline
void GNurbsSrf::GetDomainV(double &min, double &max) const
{
	min = V[q];
	max = V[n + 1];
}

inline
void GNurbsSrf::SetPrecision(double error)
{
	Precision = error;
}

inline
double GNurbsSrf::GetPrecision()
{
	return Precision;
}

