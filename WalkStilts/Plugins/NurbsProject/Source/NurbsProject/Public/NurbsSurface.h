// Copyright Zhangci 2018

#pragma once

#include "GNURBS.h"
#include "ProceduralMeshComponent.h"
#include "Components/SplineComponent.h"
#include "Components/ArrowComponent.h"
#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "NurbsSurface.generated.h"

UCLASS(BlueprintType)
class NURBSPROJECT_API ANurbsSurface : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ANurbsSurface();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

protected:

	GNurbsSrf *Surface;

	// Convert a FVector to a GVector.
	inline GVector3 FVec2GVec(FVector v) { return GVector3(v[0], v[1], v[2]); }

	// Convert a GVector to a FVector.
	inline FVector GVec2FVec(GVector3 v) { return FVector(v[0], v[1], v[2]); }

	// Convert a TArray<float> to an array of double.
	inline double* TArray2Double(TArray<float> arr)
	{
		double *tmp = new double[arr.Num()];
		for (int i = 0; i<arr.Num(); i++)
		{
			tmp[i] = arr[i];
		}
		return tmp;
	}

	// Convert a TArray<FVector> to an array of GVector3.
	inline GVector3* FVecArray2GVec3(TArray<FVector> pts)
	{
		GVector3* tmp = new GVector3[pts.Num()];
		for (int i = 0; i < pts.Num(); i++)
		{
			tmp[i] = FVec2GVec(pts[i]);
		}
		return tmp;
	}

	// Check whether the parameter is within the range of knot vector.
	inline bool CheckKnotValidity(int Degree, TArray<float> arr)
	{
		bool b = true;
		for (int i = 0; i < Degree; i++)
		{
			b &= EQ(arr[i], arr[i + 1], ERROR);
		}

		for (int i = arr.Num() - Degree - 1; i < arr.Num() - 1; i++)
		{
			b &= EQ(arr[i], arr[i + 1], ERROR);
		}

		for (int i = 0; i < arr.Num() - 1; i++)
		{
			b &= arr[i + 1] - arr[i] >= 0 ? true : false;
		}

		return b;
	}

	/*
	*  Update control points and knots of curve when curve updated.
	*/
	void Update(TArray<FVector>& Points, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV);

protected:

	UPROPERTY(VisibleAnywhere, Category = "Surface")
		UProceduralMeshComponent * SurfaceMesh;

	UPROPERTY(VisibleAnywhere, Category = "Surface")
		TArray<USplineComponent*> SplineArray;

	UPROPERTY(VisibleAnywhere, Category = "Surface")
		TArray<UArrowComponent*> NormalArray;

	/*
	*  Create a nurbs surface.
	*  DegreeU: degree for U direction.
	*  IndexU: last index of contorl points for U direction.
	*  KnotU: knot for U direction.
	*  DegreeV: degree for V direction.
	*  IndexV: last index of control pionts for V direction.
	*  KnotV: knot for V direction.
	*  Points: control points of the surface.
	*
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Nurbs Surface", Keywords = "NURBS"), Category = "NURBS|Surface")
		void CreateNurbsSurface(int DegreeU, int IndexU, TArray<float> KnotU, int DegreeV, int IndexV, TArray<float> KnotV, TArray<FVector> Points);


	/*
	*  Create a closed nurbs surface.
	*  DegreeU: degree for U direction.
	*  IndexU: last index of contorl points for U direction.
	*  DegreeV: degree for V direction.
	*  IndexV: last index of control pionts for V direction.
	*  Points: control points of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Closed Nurbs Surface", Keywords = "NURBS"), Category = "NURBS|Surface")
		void CreateNurbsClosedSurface(int DegreeU, int IndexU, int DegreeV, int IndexV, TArray<FVector> Points);

	/*
	*  Create a spline dynamiclly to draw control line of surface.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Spline", Keywords = "NURBS"), Category = "NURBS|Surface")
		USplineComponent* CreateSpline();

	/*
	*  Clear all spline lines.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Clear All Splines", Keywords = "NURBS"), Category = "NURBS|Surface")
		void ClearAllSpline();

	/*
	*  Clear all normals.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Clear All Arrows", Keywords = "NURBS"), Category = "NURBS|Surface")
		void ClearAllArrows();

	/*
	*  Get a point on the surface.
	*  u: U direction parameter, usually from 0 to 1.
	*  kth: U direction derivative order.
	*  v: V direction parameter, usually from 0 to 1.
	*  lth: V direction derivative order.
	*  Return the point on the surface according to the u and v.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Evalution Surface Point", Keywords = "NURBS"), Category = "NURBS|Surface")
		FVector Eval(float u, int kth, float v, int lth);

	/*
	*  Draw nurbs surface mesh.
	*  Vertices: vertices on surface.
	*  RowResolution: vertex count of row (u direction).
	*  ColResolution: vertex count of column (v direction).
	*  Triangles: array of the tringles' index.
	*  Normals: array of computed vertex normal.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Draw Surface Mesh", Keywords = "NURBS"), Category = "NURBS|Surface")
		void CreateMesh(int RowResolution, int ColResolution, UPARAM(ref) TArray<FVector>& Vertices, UMaterialInterface *Material, TArray<int32>& Triangles, TArray<FVector>& Normals);

	/*
	*  Insert a knot to the surface.
	*  u: u direction value to insert.
	*  r: u direction multiplicity of insertion.
	*  v: v direction value to insert.
	*  k: v direction multiplicity of insertion.
	*  knotU: u direction knot.
	*  knotV: v direction knot.
	*  Points: control points of the surface.
	*  UNum: count of u direction control points.
	*  VNum: count of v direction control points.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Insert Knots", Keywords = "NURBS"), Category = "NURBS|Surface")
		void InsertKnots(const float u, int r, const float v, int k, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV, UPARAM(ref) TArray<FVector>& Points, int& UNum, int& VNum);

	/*
	*  Remove a knot from the surface.
	*  u: u direction value to remove.
	*  r: u direction multiplicity of removal.
	*  v: v direction value to remove.
	*  k: v direction multiplicity of removal.
	*  knotU: u direction knot.
	*  knotV: v direction knot.
	*  Points: control points of the surface.
	*  UNum: count of u direction control points.
	*  VNum: count of v direction control points.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Remove Knots", Keywords = "NURBS"), Category = "NURBS|Surface")
		void RemoveKnots(const float u, int r, const float v, int k, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV, UPARAM(ref) TArray<FVector>& Points, int& UNum, int& VNum);

	/*
	*  Knots refinement.
	*  r: multiplicity of refinement.
	*  knotU: u direction knot.
	*  knotV: v direction knot.
	*  Points: control points of the surface.
	*  UNum: count of u direction control points.
	*  VNum: count of v direction control points.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Refinement", Keywords = "NURBS"), Category = "NURBS|Surface")
		void Refinement(int r, UPARAM(ref)TArray<float>& KnotU, UPARAM(ref)TArray<float>& KnotV, UPARAM(ref)TArray<FVector>& Points, int & UNum, int & VNum);


	/*
	*  Make compact form of the surface.
	*  Points: control points of the curve.
	*  knotU: u direction knot.
	*  knotV: v direction knot.
	*  UNum: count of u direction control points.
	*  VNum: count of v direction control points.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Make Compact Form", Keywords = "NURBS"), Category = "NURBS|Surface")
		void MakeCompactForm(UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV, int & UNum, int & VNum);

	/*
	*  Edit a point of the surface by giving an offset.
	*  u: parameter on u direction.
	*  v: parameter on v direction.
	*  Offset: offset vector.
	*  Points: control points of the curve.
	*  knotU: u direction knot.
	*  knotV: v direction knot.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Edit", Keywords = "NURBS"), Category = "NURBS|Surface")
		void Edit(float u, float v, FVector Offset, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV);

	/*
	*  Draw an arrow.
	*  Location: arrow's location.
	*  Direction: where is arrow pointing to.
	*  Scale: arrow's scale.
	*  Color: arrow's color.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Arrow", Keywords = "NURBS"), Category = "NURBS|Surface")
		void CreateArrorw(FVector Location, FRotator Direction, FVector Scale, FLinearColor Color);
	
	
};
