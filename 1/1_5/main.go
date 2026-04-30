package main

import (
	"fmt"
	"math"
)

func makeE(n int) [][]float64 {
	E := make([][]float64, n)
	for i := 0; i < n; i++ {
		E[i] = make([]float64, n)
		E[i][i] = 1.0
	}
	return E
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

func sign(x float64) float64 {
	if x < 0 {
		return -1.0
	}
	return 1.0
}

func copyMatrix(A [][]float64) [][]float64 {
	n := len(A)
	C := make([][]float64, n)
	for i := 0; i < n; i++ {
		C[i] = make([]float64, len(A[i]))
		copy(C[i], A[i])
	}
	return C
}

func qrDecomposition(A [][]float64) ([][]float64, [][]float64) {
	n := len(A)
	R := copyMatrix(A)
	Q := makeE(n)

	for k := 0; k < n-1; k++ {
		// v = b + sign(b_1)*||b||*e_1
		// ||b|| = sum(b_i^2)^1/2
		sum := 0.0
		for i := k; i < n; i++ {
			sum += R[i][k] * R[i][k]
		}
		normB := math.Sqrt(sum)

		v := make([]float64, n)
		for i := k; i < n; i++ {
			v[i] = R[i][k]
		}
		v[k] = v[k] + sign(R[k][k])*normB

		vTv := 0.0
		for i := k; i < n; i++ {
			vTv += v[i] * v[i]
		}

		// H = E - 2*(v*v^T)/(v^T*v)
		H := makeE(n)
		if vTv > 1e-12 {
			for i := k; i < n; i++ {
				for j := k; j < n; j++ {
					H[i][j] -= 2.0 * v[i] * v[j] / vTv
				}
			}
		}

		R = multiply(H, R)
		Q = multiply(Q, H)
	}
	return Q, R
}


func check(A [][]float64, eps float64) bool {
	n := len(A)
	for m := 0; m < n-1; m++ {
		sum := 0.0
		for l := m + 1; l < n; l++ {
			sum += A[l][m] * A[l][m]
		}

		if math.Sqrt(sum) > eps {
			return false
		}
	}
	return true
}


func getSZ(A [][]float64, eps float64) []complex128 {
	n := len(A)
	res := make([]complex128, 0, n)

	i := 0
	for i < n {
		if i == n-1 || math.Abs(A[i+1][i]) < eps {
			res = append(res, complex(A[i][i], 0))
			i++
		} else {
			a := A[i][i]
			b := A[i][i+1]
			c := A[i+1][i]
			d := A[i+1][i+1]

			tr := a + d
			det := a*d - b*c
			D := tr*tr - 4*det

			if D >= 0 {
				l1 := (tr + math.Sqrt(D)) / 2
				l2 := (tr - math.Sqrt(D)) / 2
				res = append(res, complex(l1, 0), complex(l2, 0))
			} else {
				re := tr / 2
				im := math.Sqrt(-D) / 2
				res = append(res, complex(re, im), complex(re, -im))
			}
			i += 2
		}
	}
	return res
}

func main() {
	A_k := [][]float64{
		{5, -5, -6},
		{-1, -8, -5},
		{2, 7, -3},
	}

	eps := 0.01
	n := len(A_k)
	maxIters := 50

	fmt.Printf("Start (eps = %v)\n", eps)

	k := 0
	for {
		Q_k, R_k := qrDecomposition(A_k)

		A_k = multiply(R_k, Q_k)
		k++

		if check(A_k, eps) || k >= maxIters {
			break
		}
	}

	fmt.Printf("Iterations num: %d\n\n", k)

	fmt.Println("A^k:")
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Printf("%.4f ", A_k[i][j])
		}
		fmt.Println()
	}

	sz := getSZ(A_k, eps)

	fmt.Println("\nSZ:")
	for i, lambda := range sz {
		if imag(lambda) == 0 {
			fmt.Printf("λ_%d = %.4f\n", i+1, real(lambda))
		} else {
			fmt.Printf("λ_%d = %.4f %+.4fi\n", i+1, real(lambda), imag(lambda))
		}
	}
}