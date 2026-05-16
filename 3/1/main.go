// 3.1. Используя таблицу значений   функции  , 
// вычисленных в точках     построить интерполяционные многочлены Лагранжа и Ньютона, 
// проходящие через точки  .  Вычислить значение погрешности интерполяции в точке  .

package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return math.Tan(x)
}

func lagrange(x_star float64, X []float64, Y []float64) float64 {
	n := len(X)
	var L_n float64 = 0.0

	for i := 0; i < n; i++ {
		var l_i float64 = 1.0
		for j := 0; j < n; j++ {
			if i != j {
				l_i = l_i * (x_star - X[j]) / (X[i] - X[j])
			}
		}
		L_n = L_n + Y[i]*l_i
	}
	return L_n
}

func newton(x_star float64, X []float64, Y []float64) float64 {
	n := len(X)
	
	F := make([][]float64, n)
	for i := range F {
		F[i] = make([]float64, n)
		F[i][0] = Y[i]
	}

	for j := 1; j < n; j++ {
		for i := 0; i < n-j; i++ {
			F[i][j] = (F[i+1][j-1] - F[i][j-1]) / (X[i+j] - X[i])
		}
	}

	var P_n float64 = F[0][0]
	var multiplier float64 = 1.0

	for j := 1; j < n; j++ {
		multiplier = multiplier * (x_star - X[j-1])
		P_n = P_n + F[0][j]*multiplier
	}

	return P_n
}

func printResults(X []float64, x_star float64) {
	Y := make([]float64, len(X))
	for i := 0; i < len(X); i++ {
		Y[i] = f(X[i])
	}

	fmt.Println("X_i    |     Y_i")
	for i := 0; i < len(X); i++ {
		fmt.Printf("%.5f | %.5f\n", X[i], Y[i])
	}
	fmt.Println()

	lagrange_val := lagrange(x_star, X, Y)
	newton_val := newton(x_star, X, Y)
	exact_val := f(x_star)

	fmt.Printf("X* = 3pi/16 = %.5f\n\n", x_star)
	
	fmt.Printf("L(X*) = %.10f\n", lagrange_val)
	fmt.Printf("P(X*) = %.10f\n", newton_val)
	fmt.Printf("y(X*) = %.10f\n\n", exact_val)

	fmt.Printf("Погрешность Лагранжа: %.10f\n", math.Abs(lagrange_val-exact_val))
	fmt.Printf("Погрешность Ньютона: %.10f\n", math.Abs(newton_val-exact_val))
}


func main() {
	x_star := 3.0 * math.Pi / 16.0
	X_a := []float64{
		0,
		math.Pi / 8.0,
		2.0 * math.Pi / 8.0,
		3.0 * math.Pi / 8.0,
	}
	printResults(X_a, x_star)

	X_b := []float64{
		0,
		math.Pi / 8.0,
		math.Pi / 3.0,
		3.0 * math.Pi / 8.0,
	}
	printResults(X_b, x_star)
}