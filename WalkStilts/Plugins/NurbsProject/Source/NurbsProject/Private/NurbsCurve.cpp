// Copyright Zhangci 2018

#include "NurbsCurve.h"


// Sets default values
ANurbsCurve::ANurbsCurve()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;


	bShowControlPoints = true;
	bShowPointsCoordinate = true;
	bShowKnot = true;
	bShowCurve = true;
}

// Called when the game starts or when spawned
void ANurbsCurve::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void ANurbsCurve::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}


void ANurbsCurve::Update(TArray<FVector>& Points, UPARAM(ref)TArray<float>& Knot)
{
	Points.Empty();
	Points.SetNum(Curve->GetCtlPtIdx() + 1);

	for (int i = 0; i < Points.Num(); i++)
	{
		Points[i] = GVec2FVec(Curve->GetCtlPt()[i]);
	}

	Knot.Empty();
	Knot.SetNum(Curve->GetKnots().GetKnotIdx() + 1);

	for (int i = 0; i < Knot.Num(); i++)
	{
		Knot[i] = (float)Curve->GetKnots()[i];
	}
}

void ANurbsCurve::CreateText(FVector Position, FText Text, float Size)
{
	UTextRenderComponent* text = NewObject<UTextRenderComponent>(this);

	text->SetText(Text);
	text->SetRelativeLocation(Position);
	text->SetRelativeRotation(FRotator(0, 180, 0));
	text->SetHorizontalAlignment(EHorizTextAligment::EHTA_Center);
	text->SetWorldSize(Size);
	text->SetupAttachment(RootComponent);
	text->RegisterComponent();

}

void ANurbsCurve::CreateNurbsCurve(int Degree, TArray<FVector> Points, TArray<float> Knot)
{
	if (Knot.Num() - 1 != Points.Num() + Degree)
	{
		UE_LOG(LogClass, Error, TEXT("CreateNurbsCurve: m must be n + P + 1."));
		return;
	}

	if (CheckKnotValidity(Degree, Knot))
		Curve = new GNurbsCrv(Degree, Points.Num() - 1, FVecArray2GVec3(Points), TArray2Double(Knot));
	else
	{
		UE_LOG(LogClass, Error, TEXT("CreateNurbsCurve: the knot vector must be {a, ..., a, U_p+1, ..., U_m-p-1, b, ..., b}, and the count of a,b is Degree + 1."));
		return;
	}
}

void ANurbsCurve::CreateNurbsClosedCurve(int Degree, TArray<FVector> Points)
{
	Curve = get_gnurbs_closed_crv(Degree, Points.Num() - 1, FVecArray2GVec3(Points));
}

FVector ANurbsCurve::Eval(float u, int nth)
{
	if (Curve == nullptr)
		return FVector(0, 0, 0);

	return GVec2FVec(Curve->Eval(u, nth));
}

void ANurbsCurve::InsertKnots(const float u, int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot)
{

	if (Curve == nullptr)
		return;

	Curve->InsertKnots(u, r);
	Update(Points, Knot);
}

int ANurbsCurve::RemoveKnots(const float u, int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot)
{
	if (Curve == nullptr)
		return -1;

	int tmp = Curve->RemoveKnots(u, r);
	Update(Points, Knot);
	return tmp;
}

void ANurbsCurve::Refinement(const int r, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot)
{
	if (Curve == nullptr)
		return;

	for (int i = 0; i < r; i++)
		Curve->Refinement();

	Update(Points, Knot);
}

void ANurbsCurve::Edit(float u, FVector offset, UPARAM(ref) TArray<FVector>& Points, UPARAM(ref) TArray<float>& Knot)
{
	if (Curve == nullptr)
		return;

	Curve->Edit(u, FVec2GVec(offset));

	Update(Points, Knot);

}

void ANurbsCurve::MakeBezierForm(UPARAM(ref)TArray<FVector>& Points, UPARAM(ref)TArray<float>& Knot)
{
	if (Curve == nullptr)
		return;

	Curve->MakeBezierForm();
	Update(Points, Knot);
}

void ANurbsCurve::MakeCompactForm(UPARAM(ref)TArray<FVector>& Points, UPARAM(ref)TArray<float>& Knot)
{
	if (Curve == nullptr)
		return;

	Curve->MakeCompactForm();
	Update(Points, Knot);
}

