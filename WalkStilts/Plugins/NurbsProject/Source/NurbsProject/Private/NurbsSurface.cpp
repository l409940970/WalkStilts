// Copyright Zhangci 2018
#include "NurbsSurface.h"


// Sets default values
ANurbsSurface::ANurbsSurface()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	SurfaceMesh = CreateDefaultSubobject<UProceduralMeshComponent>(TEXT("Surface Mesh"));

	SurfaceMesh->AttachToComponent(RootComponent, FAttachmentTransformRules::KeepRelativeTransform);
}

// Called when the game starts or when spawned
void ANurbsSurface::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void ANurbsSurface::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

void ANurbsSurface::Update(TArray<FVector>& Points, UPARAM(ref)TArray<float>& KnotU, UPARAM(ref)TArray<float>& KnotV)
{
	Points.Empty();

	for (int i = 0; i <= Surface->GetCtlPtIdxU(); i++)
	{
		for (int j = 0; j <= Surface->GetCtlPtIdxV(); j++)
		{
			Points.Add(GVec2FVec(Surface->GetCtlPt()[i][j]));
		}
	}

	KnotU.Empty();
	KnotU.SetNum(Surface->GetKnotU().GetKnotIdx() + 1);

	KnotV.Empty();
	KnotV.SetNum(Surface->GetKnotV().GetKnotIdx() + 1);


	for (int i = 0; i < KnotU.Num(); i++)
	{
		KnotU[i] = (float)Surface->GetKnotU()[i];
	}

	for (int i = 0; i < KnotV.Num(); i++)
	{
		KnotV[i] = (float)Surface->GetKnotU()[i];
	}
}

void ANurbsSurface::CreateNurbsSurface(int DegreeU, int IndexU, TArray<float> KnotU, int DegreeV, int IndexV, TArray<float> KnotV, TArray<FVector> Points)
{
	if (KnotU.Num() - 1 != IndexU + DegreeU + 1)
	{
		UE_LOG(LogClass, Error, TEXT("CreateNurbsSurface: m must be n + P + 1.(U direction)"));
		return;
	}

	if (KnotV.Num() - 1 != IndexV + DegreeV + 1)
	{
		UE_LOG(LogClass, Error, TEXT("CreateNurbsSurface: m must be n + P + 1.(V direction)"));
		return;
	}

	if (CheckKnotValidity(DegreeU, KnotU) && CheckKnotValidity(DegreeV, KnotV))
		Surface = new GNurbsSrf(DegreeU, IndexU, TArray2Double(KnotU), DegreeV, IndexV, TArray2Double(KnotV), FVecArray2GVec3(Points));
	else
	{
		UE_LOG(LogClass, Error, TEXT("CreateNurbsSurface: the knot vector must be {a, ..., a, U_p+1, ..., U_m-p-1, b, ..., b}, and the count of a,b is Degree + 1."));
		return;
	}


}

void ANurbsSurface::CreateNurbsClosedSurface(int DegreeU, int IndexU, int DegreeV, int IndexV, TArray<FVector> Points)
{
	Surface = get_gnurbs_closed_srf(DegreeU, IndexU, DegreeV, IndexV, FVecArray2GVec3(Points), true, true);
}

USplineComponent * ANurbsSurface::CreateSpline()
{

	USplineComponent *spline = NewObject<USplineComponent>(this);
	spline->RegisterComponent();
	spline->AttachToComponent(RootComponent, FAttachmentTransformRules::KeepRelativeTransform);
	spline->ClearSplinePoints();
	spline->SetUnselectedSplineSegmentColor(FLinearColor(1, 0, 0, 1));

	SplineArray.Add(spline);

	return spline;

}

void ANurbsSurface::ClearAllSpline()
{
	for (int i = 0; i < SplineArray.Num(); ++i)
	{
		if (SplineArray[i] != nullptr)
		{
			SplineArray[i]->ClearSplinePoints();
			SplineArray[i]->DestroyComponent();
		}
	}

	SplineArray.Empty();
}

void ANurbsSurface::ClearAllArrows()
{
	for (int i = 0; i < NormalArray.Num(); ++i)
	{
		if (NormalArray[i] != nullptr)
		{
			NormalArray[i]->DestroyComponent();
		}
	}

	NormalArray.Empty();
}

FVector ANurbsSurface::Eval(float u, int kth, float v, int lth)
{
	if (Surface == nullptr)
		return FVector(0, 0, 0);

	return GVec2FVec(Surface->Eval(u, v, kth, lth));
}

void ANurbsSurface::CreateMesh(int RowResolution, int ColResolution, UPARAM(ref) TArray<FVector>& Vertices, UMaterialInterface *Material, TArray<int32>& Triangles, TArray<FVector>& Normals)
{

	/* Triangles
	*
	*  Compute index of all triangles.
	*   3---2
	*   |   |
	*   0---1
	*   the face triangles [0 2 1] and [0 3 2] are computed.
	*/

	TArray<int32> triangles;

	for (int i = 0; i < RowResolution - 1; i++)
	{
		for (int j = 0; j < ColResolution - 1; j++)
		{
			int v0 = i*RowResolution + j;
			int v1 = i*RowResolution + j + 1;
			int v2 = (i + 1)*ColResolution + (j + 1);
			int v3 = (i + 1)*ColResolution + j;

			triangles.Add(v0);
			triangles.Add(v2);
			triangles.Add(v1);

			triangles.Add(v0);
			triangles.Add(v3);
			triangles.Add(v2);
		}
	}

	Triangles = triangles;

	/* Normals
	*
	*  Compute normals of all vertices.
	*       3
	*       |
	*   1---0---2
	*       |
	*       4
	*  So, normal of vertex 0 will be: (1-0)x(4-0) + (4-0)x(2-0) + (2-0)x(3-0) + (3-0)x(1-0)
	*  Note: this method is one of the way to compute vertices' normal. And you also can compute normals when computing triangles.
	*/

	TArray<FVector> normals;

	for (int i = 0; i < RowResolution; i++)
	{
		for (int j = 0; j < ColResolution; j++)
		{
			FVector L(0, 0, 0), R(0, 0, 0), U(0, 0, 0), D(0, 0, 0);
			//int nL = 

			if (j - 1 >= 0)
				L = Vertices[i*RowResolution + j - 1] - Vertices[i*RowResolution + j];
			if (j + 1 < ColResolution)
				R = Vertices[i*RowResolution + j + 1] - Vertices[i*RowResolution + j];
			if (i + 1 < RowResolution)
				U = Vertices[(i + 1)*ColResolution + j] - Vertices[i*RowResolution + j];
			if (i - 1 >= 0)
				D = Vertices[(i - 1)*ColResolution + j] - Vertices[i*RowResolution + j];

			FVector v = FVector::CrossProduct(L, D) + FVector::CrossProduct(D, R) + FVector::CrossProduct(R, U) + FVector::CrossProduct(U, L);
			normals.Add(v.GetSafeNormal());

		}
	}

	Normals = normals;

	/*
	*  Texture coordinate.
	*/

	TArray<FVector2D> UV0;
	for (int i = 0; i < RowResolution; i++)
	{
		for (int j = 0; j < ColResolution; j++)
		{
			UV0.Add(FVector2D((float)j / (ColResolution - 1), (float)i / (RowResolution - 1)));
		}
	}

	/*
	*  Vertex color.
	*/

	TArray<FColor> vertexColors;
	for (int i = 0; i < Vertices.Num(); i++)
	{
		vertexColors.Add(FColor(255, 255, 255));
	}

	/*
	*  Tangents.
	*/
	TArray<FProcMeshTangent> tangents;
	for (int i = 0; i < Vertices.Num(); i++)
	{
		tangents.Add(FProcMeshTangent(0, 1, 0));
	}

	SurfaceMesh->ClearAllMeshSections();
	SurfaceMesh->CreateMeshSection(1, Vertices, triangles, normals, UV0, vertexColors, tangents, true);
	SurfaceMesh->SetMaterial(1, Material);
	SurfaceMesh->SetCollisionEnabled(ECollisionEnabled::QueryAndPhysics);

}

void ANurbsSurface::InsertKnots(const float u, int r, const float v, int k, UPARAM(ref)TArray<float>& KnotU, UPARAM(ref)TArray<float>& KnotV, UPARAM(ref)TArray<FVector>& Points, int& UNum, int& VNum)
{
	if (Surface == nullptr)
		return;

	Surface->InsertKnotsU(u, r);
	Surface->InsertKnotsV(v, k);

	Update(Points, KnotU, KnotV);

	UNum = Surface->GetCtlPtIdxU() + 1;
	VNum = Surface->GetCtlPtIdxV() + 1;
}

void ANurbsSurface::RemoveKnots(const float u, int r, const float v, int k, UPARAM(ref)TArray<float>& KnotU, UPARAM(ref)TArray<float>& KnotV, UPARAM(ref)TArray<FVector>& Points, int & UNum, int & VNum)
{
	if (Surface == nullptr)
		return;

	Surface->RemoveKnotsU(u, r);
	Surface->RemoveKnotsV(v, k);

	Update(Points, KnotU, KnotV);

	UNum = Surface->GetCtlPtIdxU() + 1;
	VNum = Surface->GetCtlPtIdxV() + 1;
}

void ANurbsSurface::Refinement(int r, UPARAM(ref)TArray<float>& KnotU, UPARAM(ref)TArray<float>& KnotV, UPARAM(ref)TArray<FVector>& Points, int & UNum, int & VNum)
{
	if (Surface == nullptr)
		return;

	for (int i = 0; i < r; i++)
	{
		Surface->RefinementU();
		Surface->RefinementV();
	}

	Update(Points, KnotU, KnotV);

	UNum = Surface->GetCtlPtIdxU() + 1;
	VNum = Surface->GetCtlPtIdxV() + 1;
}

void ANurbsSurface::MakeCompactForm(UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV, int & UNum, int & VNum)
{
	if (Surface == nullptr)
		return;

	Surface->MakeCompactForm();

	Update(Points, KnotU, KnotV);

	UNum = Surface->GetCtlPtIdxU() + 1;
	VNum = Surface->GetCtlPtIdxV() + 1;
}

void ANurbsSurface::Edit(float u, float v, FVector Offset, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& KnotU, UPARAM(ref) TArray<float>& KnotV)
{
	if (Surface == nullptr)
		return;

	Surface->Edit(u, v, FVec2GVec(Offset));

	Update(Points, KnotU, KnotV);

}

void ANurbsSurface::CreateArrorw(FVector Location, FRotator Direction, FVector Scale, FLinearColor Color)
{

	UArrowComponent* arrow = NewObject<UArrowComponent>(this);
	arrow->RegisterComponent();
	arrow->AttachToComponent(RootComponent, FAttachmentTransformRules::KeepRelativeTransform);
	arrow->SetHiddenInGame(false);
	arrow->SetArrowColor(Color);
	arrow->SetRelativeTransform(FTransform(Direction, Location, Scale));

	NormalArray.Add(arrow);

}
