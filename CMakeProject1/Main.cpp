#include <iostream>
#include <vector>
#include <numeric>
#include <cassert>

typedef uint32_t uint32;
typedef std::vector<float> Vector;

const float kTolerance = 1.0e-10f;
const float kToleranceSquared = kTolerance * kTolerance;

const float A[] = { 1.0f, -2.0f, 1.0f };


void Print(const Vector& V)
{
	for (float c : V)
	{
		printf("%.f ", c);
	}

	printf("\n");
}


Vector Add(const Vector& U, const Vector& V)
{
	Vector result(U.size());
	for (uint32 i = 0; i < U.size(); i++)
	{
		result[i] = U[i] + V[i];
	}
	return result;
}

Vector Multiply(const Vector& U, const Vector& V)
{
	Vector result(U.size());
	for (uint32 i = 0; i < U.size(); i++)
	{
		result[i] = U[i] * V[i];
	}
	return result;
}

Vector Multiply(const Vector& V, float value)
{
	Vector result(V.size());
	for (uint32 i = 0; i < V.size(); i++)
	{
		result[i] = V[i] * value;
	}
	return result;
}

float DotProduct(const Vector& U, const Vector& V)
{
	float result = 0.0f;
	for (uint32 i = 0; i < U.size(); i++)
	{
		result += U[i] * V[i];
	}
	return result;
}

float MaxComponent(const Vector& U)
{
	float max = FLT_MIN;
	for (float c : U)
	{
		max = std::max(max, c);
	}
	return max;
}

float Length2(const Vector& V)
{
	return DotProduct(V, V);
}

float Length(const Vector& V)
{
	return std::sqrt(Length2(V));
}

Vector MultiplyMatrixVector(const Vector& U)
{
	uint32 N = U.size();
	Vector Result(N, 0.0f);
	Result[0] += A[1] * U[0];
	Result[0] += A[2] * U[1];

	for (uint32 i = 0; i < N - 2; i++)
	{
		Result[i+1] += A[0] * U[i + 0];
		Result[i+1] += A[1] * U[i + 1];
		Result[i+1] += A[2] * U[i + 2];
	}

	Result[N - 1] += A[1] * U[N - 2];
	Result[N - 1] += A[2] * U[N - 1];
	return Result;
}

Vector ConjugateGradientSolver(const Vector& B)
{
	Vector X(B.size(), 0.0);
	Vector R = B;
	Vector P = R;
	Vector ROld = R;

	float RR = Length2(R);
	for (uint32 k = 0; k < B.size(); k++)
	{
		Vector AP = MultiplyMatrixVector(P);

		float RROld = RR;
		float Alpha = RR / DotProduct(P, AP);
		assert(std::isfinite(Alpha));

		X = Add(X, Multiply(P, Alpha));
		R = Add(R, Multiply(AP, -Alpha));

		if (Length(R) < kToleranceSquared)
		{
			break;
		}

		RR = Length2(R);

		float Beta = RR / RROld;
		assert(std::isfinite(Beta));

		P = Add(R, Multiply(P, Beta));
	}

	return X;
}

int main(int argc, char* argv[])
{
	// TODO: Read from file
	uint32 N = 3;
	Vector B(N);
	for (uint32 i = 0; i < N; i++)
	{
		B[i] = 0.0f;
	}

	B[(N - 1) / 2] = 1.0f;

	Vector Result = ConjugateGradientSolver(B);
	Print(Result);
}
