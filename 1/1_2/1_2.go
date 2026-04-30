package main

import "fmt"

// b - главная диагональ
func Algorithm(a, b, c, d []float64) []float64 {
	n := len(d)
	P := make([]float64, n)
	Q := make([]float64, n)
	x := make([]float64, n)

	P[0] = -c[0] / b[0]
	Q[0] = d[0] / b[0]

	for i := 1; i < n; i++ {
		denominator := a[i]*P[i-1] + b[i]

		if i < n-1 {
			P[i] = -c[i] / denominator
		}
		
		Q[i] = (d[i] - a[i]*Q[i-1]) / denominator
	}

	x[n-1] = Q[n-1]
	
	for i := n - 2; i >= 0; i-- {
		x[i] = P[i]*x[i+1] + Q[i]
	}

	return x
}

func main() {
	a := []float64{0, -4, -1, 6, 4}
	b := []float64{13, 9, -12, 20, 5}
	c := []float64{-5, -5, -6, -5, 0}
	d := []float64{-66, -47, -43, -74, 14} 

	solution := Algorithm(a, b, c, d)

	for i, val := range solution {
		fmt.Printf("x_%d = %.1f\n", i+1, val)
	}
}
