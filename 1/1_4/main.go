package main

import (
	"fmt"
	"math"
)

const EPS = 0.001

func calcT(A [][]float64) float64 {
	sum := 0.0
	n := len(A)
	for l := 0; l < n; l++ {
		for m := l + 1; m < n; m++ {
			sum += A[l][m] * A[l][m]
		}
	}
	return math.Sqrt(sum)
}

func findMaxOffDiagonal(A [][]float64) (float64, int, int) {
	n := len(A)
	var maxVal float64 = -1.0
	var maxI, maxJ int

	for l := 0; l < n; l++ {
		for m := l + 1; m < n; m++ {
			val := math.Abs(A[l][m])
			if val > maxVal {
				maxVal = val
				maxI = l
				maxJ = m
			}
		}
	}
	return A[maxI][maxJ], maxI, maxJ
}


func multiply(A, B [][]float64) [][]float64 {
	n := len(A)
	C := make([][]float64, n)
	for i := 0; i < n; i++ {
		C[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}
	return C
}

func transpose(A [][]float64) [][]float64 {
	n := len(A)
	T := make([][]float64, n)
	for i := 0; i < n; i++ {
		T[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			T[i][j] = A[j][i]
		}
	}
	return T
}

func makeE(n int) [][]float64 {
	E := make([][]float64, n)
	for i := 0; i < n; i++ {
		E[i] = make([]float64, n)
		E[i][i] = 1.0
	}
	return E
}

func main() {
	A := [][]float64{
		{5, 5, 3},
		{5, -4, 1},
		{3, 1, 2},
	}

	n := len(A)
	U := makeE(n)
	
	k := 0
	for {
		tA := calcT(A)
		fmt.Printf("iter %d: t(A) = %.4f\n", k, tA)

		if tA < EPS {
			break
		}

		a_ij, i, j := findMaxOffDiagonal(A)

		var phi float64
		if math.Abs(A[i][i] - A[j][j]) < 1e-12 {
			phi = math.Pi / 4.0
		} else {
			phi = 0.5 * math.Atan((2.0 * a_ij)/(A[i][i] - A[j][j]))
		}

		U_k := makeE(n)
		U_k[i][i] = math.Cos(phi)
		U_k[j][j] = math.Cos(phi)
		U_k[i][j] = -math.Sin(phi)
		U_k[j][i] = math.Sin(phi)

		U_k_T := transpose(U_k)
		temp := multiply(U_k_T, A)
		A = multiply(temp, U_k)

		U = multiply(U, U_k)

		k++
	}

	fmt.Printf("Total iter num: %d\n\n", k)

	fmt.Println("СЗ:")
	for i := 0; i < n; i++ {
		fmt.Printf("λ_%d ≈ %.4f\n", i+1, A[i][i])
	}

	fmt.Println("\nU:")
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Printf("%.4f ", U[i][j])
		}
		fmt.Println()
	}
}