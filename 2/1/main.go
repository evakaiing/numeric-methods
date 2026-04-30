package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return math.Sqrt(1-x*x) - math.Exp(x) + 0.1
}

func df(x float64) float64 {
	return -x/math.Sqrt(1-x*x) - math.Exp(x)
}

func phi(x float64) float64 {
	return x + 0.9*(math.Sqrt(1-x*x)-math.Exp(x)+0.1)
}

func newtonMethod(x0, eps float64) (float64, int) {
	xk := x0
	k := 0

	for {
		// x_new = x - f(x)/f'(x)
		xNext := xk - f(xk)/df(xk)
		
		err := math.Abs(xNext - xk)
		k++
		
		if err < eps {
			return xNext, k
		}
		xk = xNext
	}
}

func simpleIterationMethod(x0, eps, q float64) (float64, int) {
	xk := x0
	k := 0
	
	factor := q / (1 - q)

	for {
		// x_new = phi(x)
		xNext := phi(xk)
		
		err := factor * math.Abs(xNext - xk)
		k++
		
		if err < eps {
			return xNext, k
		}
		xk = xNext
	}
}

func main() {
	x0_newton := 0.15
	x0_simple := 0.1
	q := 0.055

	epsilons := []float64{1e-2, 1e-4, 1e-6, 1e-8}

	for _, eps := range epsilons {
		rootNewton, itersNewton := newtonMethod(x0_newton, eps)
		rootSimple, itersSimple := simpleIterationMethod(x0_simple, eps, q)

		fmt.Printf("Точность: %v\n", eps)
		fmt.Printf("Метод Ньютона: корень = %.8f, итераций = %d\n", rootNewton, itersNewton)
		fmt.Printf("Метод простой итерации: корень = %.8f, итераций = %d\n\n", rootSimple, itersSimple)
	}
}