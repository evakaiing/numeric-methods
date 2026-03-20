package main

import (
	"fmt"
	"math"
)

const (
	maxIterCount = 1000
	epsilon      = 1e-6
)

// SeidelIterations решает СЛАУ методом Зейделя.
func SeidelIterations(A [][]float64, b []float64, epsilon float64) ([]float64, int, error) {
	n := len(A)

	for col := 0; col < n; col++ {
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
			b[col], b[pivotRow] = b[pivotRow], b[col]
		}
	}

	// x = alpha*x + beta
	alpha := make([][]float64, n)
	beta := make([]float64, n)
	for row := 0; row < n; row++ {
		alpha[row] = make([]float64, n)

		beta[row] = b[row] / A[row][row]
		for col := 0; col < n; col++ {
			if row != col {
				alpha[row][col] = -A[row][col] / A[row][row]
			}
		}
	}

	// c
	normAlpha := 0.0
	for row := 0; row < n; row++ {
		rowSum := 0.0
		for col := 0; col < n; col++ {
			rowSum += math.Abs(alpha[row][col])
		}
		normAlpha = math.Max(normAlpha, rowSum)
	}

	fmt.Printf("||alpha||_c = %.4f\n", normAlpha)

	// ||alpha|| / (1 - ||alpha||)
	coef := 1.0
	if normAlpha < 1 {
		coef = normAlpha / (1 - normAlpha)
	}

	x := copyVector(beta)
	xOld := make([]float64, n)

	for k := 1; k <= maxIterCount; k++ {
		copy(xOld, x)

		for row := 0; row < n; row++ {
			sum := beta[row]
			for col := 0; col < n; col++ {
				// if col >= row - old
				// else - new
				sum += alpha[row][col] * x[col]
			}
			x[row] = sum // сразу
		}

		// ||x(k) - x(k-1)||_c
		diffNorm := 0.0
		for row := 0; row < n; row++ {
			diffNorm = math.Max(diffNorm, math.Abs(x[row]-xOld[row]))
		}

		if coef*diffNorm <= epsilon {
			return x, k, nil
		}
	}

	return x, maxIterCount, fmt.Errorf("max iterations reached")
}

func copyVector(b []float64) []float64 {
	cp := make([]float64, len(b))
	copy(cp, b)
	return cp
}

func main() {
	A := [][]float64{
		{-23, -7, 5, 2},
		{-7, -21, 4, 9},
		{9, 5, -31, -8},
		{0, 1, -2, 10},
	}
	b := []float64{-26, -55, -58, -24}

	fmt.Println("\nМЕТОД ЗЕЙДЕЛЯ")
	x, iters, err := SeidelIterations(A, b, epsilon)
	if err != nil {
		fmt.Printf("Error: %v\n", err)
	}

	fmt.Printf("%d iterations:\n", iters)
	for row, val := range x {
		fmt.Printf("x_%d = %.4f\n", row+1, val)
	}
}
