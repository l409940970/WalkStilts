// Copyright Zhangci 2018

#pragma once

#include "GNURBS.h"
#include "GameFramework/Actor.h"
#include "Components/TextRenderComponent.h"
#include "Components/SplineComponent.h"
#include "CoreMinimal.h"
#include "NurbsCurve.generated.h"

UCLASS()
class NURBSPROJECT_API ANurbsCurve : public AActor
{
	GENERATED_BODY()
	
public:	
	// Sets default values for this actor's properties
	ANurbsCurve();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

private:

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
	void Update(TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

protected:
	GNurbsCrv *Curve;

	UPROPERTY(EditAnywhere, meta = (DisplayName = "Show Control Points", Keywords = "Coordinate"), Category = "NURBS")
		bool bShowControlPoints;

	UPROPERTY(EditAnywhere, meta = (DisplayName = "Show Points Coordinate", Keywords = "Coordinate"), Category = "NURBS")
		bool bShowPointsCoordinate;

	UPROPERTY(EditAnywhere, meta = (DisplayName = "Show Knot Vector", Keywords = "Knot"), Category = "NURBS")
		bool bShowKnot;

	UPROPERTY(EditAnywhere, meta = (DisplayName = "Show Curve", Keywords = "Curve"), Category = "NURBS")
		bool bShowCurve;

protected:

	/*
	* Create text to show.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Text", Keywords = "NURBS"), Category = "NURBS")
		void CreateText(FVector Position, FText Text, float Size = 20.0f);

	/*
	*  Create a nurbs curve.
	*  Degree: degree of the curve.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*  Note: if degree = p, last index of control points is n, last index of knot is m, then m = n + p + 1.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Nurbs Curve", Keywords = "NURBS"), Category = "NURBS|Curve")
		void CreateNurbsCurve(int Degree, TArray<FVector> Points, TArray<float> Knot);


	/*
	*  Create a closed nurbs curve.
	*  Degree: degree of the curve.
	*  Points: control points of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Create Closed Nurbs Curve", Keywords = "NURBS"), Category = "NURBS|Curve")
		void CreateNurbsClosedCurve(int Degree, TArray<FVector> Points);

	/*
	*  Get a point on the curve.
	*  u: parameter, usually from 0 to 1.
	*  nth: derivative order.
	*  Return the point on the curve according to the u.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Evalution Curve Point", Keywords = "NURBS"), Category = "NURBS|Curve")
		FVector Eval(float u, int nth);

	/*
	*  Insert a knot to the curve.
	*  u: value to insert.
	*  r: multiplicity of insertion.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Insert Knots", Keywords = "NURBS"), Category = "NURBS|Curve")
		void InsertKnots(const float u, int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

	/*
	*  Remove the knot value.
	*  u: knot value to remove.
	*  r: multiplicity of removal.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*  Return the actual multiplicity of removal.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Remove Knots", Keywords = "NURBS"), Category = "NURBS|Curve")
		int	RemoveKnots(const float u, int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

	/*
	*  Knots refinement.
	*  r: multiplicity of refinement.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Refinement", Keywords = "NURBS"), Category = "NURBS|Curve")
		void Refinement(const int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

	/*
	*  Edit a point of the curve by giving an offset.
	*  u: parameter.
	*  offset: offset vector.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Edit", Keywords = "NURBS"), Category = "NURBS|Curve")
		void Edit(float u, FVector offset, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

	/*
	*  Make bezier form of the curve.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Make Bezier Form", Keywords = "NURBS"), Category = "NURBS|Curve")
		void MakeBezierForm(UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);

	/*
	*  Make compact form of the curve.
	*  Points: control points of the curve.
	*  Knot: knot vector of the curve.
	*/
	UFUNCTION(BlueprintCallable, meta = (DisplayName = "Make Compact Form", Keywords = "NURBS"), Category = "NURBS|Curve")
		void MakeCompactForm(UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot);



	// Getters

	UFUNCTION(BlueprintPure, meta = (CompactNodeTitle = "Show Points", Keywords = "NURBS"), Category = "NURBS|Curve")
		FORCEINLINE bool GetShowPoints() { return bShowControlPoints; }

	UFUNCTION(BlueprintPure, meta = (CompactNodeTitle = "Show Coordinate", Keywords = "NURBS"), Category = "NURBS|Curve")
		FORCEINLINE bool GetShowCoordinate() { return bShowPointsCoordinate; }

	UFUNCTION(BlueprintPure, meta = (CompactNodeTitle = "Show Knot", Keywords = "NURBS"), Category = "NURBS|Curve")
		FORCEINLINE bool GetShowKnot() { return bShowKnot; }

	UFUNCTION(BlueprintPure, meta = (CompactNodeTitle = "Show Curve", Keywords = "NURBS"), Category = "NURBS|Curve")
		FORCEINLINE bool GetShowCurve() { return bShowCurve; }
	
	
};
