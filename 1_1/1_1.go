package main

import (
	"fmt"
	"math"
)

func LUDecomposition(A [][]float64) ([]int, int, error) {
	n := len(A)
	P := make([]int, n)
	for col := range P {
		P[col] = col
	}
	swaps := 0

	for col := 0; col < n; col++ {
		// диагональное преобладание
		pivotRow := col
		maxVal := math.Abs(A[col][col])
		for row := col + 1; row < n; row++ {
			if math.Abs(A[row][col]) > maxVal {
				maxVal = math.Abs(A[row][col])
				pivotRow = row
			}
		}

		if maxVal < 1e-12 {
			return nil, 0, fmt.Errorf("det == 0")
		}

		if pivotRow != col {
			A[col], A[pivotRow] = A[pivotRow], A[col]
			P[col], P[pivotRow] = P[pivotRow], P[col]
			swaps++
		}

		// прямой ход
		for row := col + 1; row < n; row++ {
			A[row][col] /= A[col][col] // коэффициент l_{row,col}

			for k := col + 1; k < n; k++ { // по столбцам правее
				A[row][k] -= A[row][col] * A[col][k]
			}
		}
	}

	return P, swaps, nil
}

func LUSolve(LU [][]float64, P []int, b []float64) []float64 {
	n := len(LU)
	x := make([]float64, n)
	z := make([]float64, n)

	for col := 0; col < n; col++ {
		z[col] = b[P[col]]
	}

	// L z = b
	for row := 0; row < n; row++ {
        for col := 0; col < row; col++ {
            z[row] -= LU[row][col] * z[col]
        }
    }

	// U x = z
	for row := n - 1; row >= 0; row-- {
        x[row] = z[row]
        for col := row + 1; col < n; col++ {
            x[row] -= LU[row][col] * x[col]
        }
        x[row] /= LU[row][row]
    }


	return x
}

func Determinant(LU [][]float64, swaps int) float64 {
	det := 1.0
	for col := 0; col < len(LU); col++ {
		det *= LU[col][col]
	}
	if swaps % 2 != 0 {
		det = -det
	}
	return det
}

func InverseMatrix(LU [][]float64, P []int) [][]float64 {
	n := len(LU)
	inv := make([][]float64, n)
	for col := range inv {
		inv[col] = make([]float64, n)
	}

	for col := 0; col < n; col++ {
		e := make([]float64, n)
		e[col] = 1.0
		x := LUSolve(LU, P, e)

		for row := 0; row < n; row++ {
			inv[row][col] = x[row]
		}
	}

	return inv
}

func copzMatrix(A [][]float64) [][]float64 {
	n := len(A)
	copyA := make([][]float64, n)
	for col := range A {
		copyA[col] = make([]float64, len(A[col]))
		copy(copyA[col], A[col])
	}
	return copyA
}
func main() {
	A := [][]float64{
		{9, -5, -6, 3},
		{1, -7, 1, 0},
		{3, -4, 9, 0},
		{6, -1, 9, 8},
	}
	b := []float64{-8, 38, 47, -8}

	LU := copzMatrix(A)
	P, swaps, err := LUDecomposition(LU)
	if err != nil {
		fmt.Println(err)
		return
	}

	x := LUSolve(LU, P, b)
	fmt.Printf("x = %v\n", x)

	det := Determinant(LU, swaps)
	fmt.Printf("det(A) = %.4f\n", det)

	inv := InverseMatrix(LU, P)
	fmt.Println("A^{-1}:")
	for _, row := range inv {
		for _, val := range row {
			fmt.Printf("%.4f ", val)
		}
		fmt.Println()
	}
}
