package main

import (
	"fmt"
	"math"
)

func f1(x1, x2 float64) float64 {
	return (x1*x1+16)*x2 - 64
}

func f2(x1, x2 float64) float64 {
	return (x1-2)*(x1-2) + (x2-2)*(x2-2) - 16
}

func df1_dx1(x1, x2 float64) float64 { return 2 * x1 * x2 }
func df1_dx2(x1, x2 float64) float64 { return x1*x1 + 16 }
func df2_dx1(x1, x2 float64) float64 { return 2 * (x1 - 2) }
func df2_dx2(x1, x2 float64) float64 { return 2 * (x2 - 2) }

func det(matrix [2][2]float64) float64 {
	return matrix[0][0]*matrix[1][1] - matrix[0][1]*matrix[1][0]
}

func phi1(x2 float64) float64 {
	return 2.0 + math.Sqrt(16.0-(x2-2.0)*(x2-2.0))
}

func phi2(x1 float64) float64 {
	return 64.0 / (x1*x1 + 16.0)
}

func dphi1_dx1(x2 float64) float64 { return 0.0 }
func dphi1_dx2(x2 float64) float64 {
	return -(x2 - 2.0) / math.Sqrt(16.0-(x2-2.0)*(x2-2.0))
}
func dphi2_dx1(x1 float64) float64 {
	return -128.0 * x1 / ((x1*x1 + 16.0) * (x1*x1 + 16.0))
}
func dphi2_dx2(x1 float64) float64 { return 0.0 }

func newtonMethod(x1_0, x2_0, eps float64) (float64, float64, int, error) {
	x1 := x1_0
	x2 := x2_0
	k := 0

	for {
		F1 := f1(x1, x2)
		F2 := f2(x1, x2)

		J := [2][2]float64{
			{df1_dx1(x1, x2), df1_dx2(x1, x2)},
			{df2_dx1(x1, x2), df2_dx2(x1, x2)},
		}

		detJ := det(J)

		if math.Abs(detJ) < 1e-12 {
			return x1, x2, k, fmt.Errorf("det(J) == 0")
		}

		A1 := [2][2]float64{
			{F1, J[0][1]},
			{F2, J[1][1]},
		}
		A2 := [2][2]float64{
			{J[0][0], F1},
			{J[1][0], F2},
		}

		x1Next := x1 - det(A1)/detJ
		x2Next := x2 - det(A2)/detJ

		errCalc := math.Max(math.Abs(x1Next-x1), math.Abs(x2Next-x2))
		k++

		if errCalc < eps {
			return x1Next, x2Next, k, nil
		}

		x1 = x1Next
		x2 = x2Next
	}
}

func simpleIterationMethod(x1_0, x2_0, eps float64) (float64, float64, int, error) {
	x1 := x1_0
	x2 := x2_0
	k := 0

	for {
		norm1 := math.Abs(dphi1_dx1(x2)) + math.Abs(dphi1_dx2(x2))
		norm2 := math.Abs(dphi2_dx1(x1)) + math.Abs(dphi2_dx2(x1))
		
		q := math.Max(norm1, norm2)

		if q >= 1.0 {
			return x1, x2, k, fmt.Errorf("(q = %.3f >= 1), iter %d", q, k)
		}

		factor := q / (1.0 - q)

		x1Next := phi1(x2)
		x2Next := phi2(x1)

		errCalc := factor * math.Max(math.Abs(x1Next-x1), math.Abs(x2Next-x2))
		k++

		if errCalc < eps {
			return x1Next, x2Next, k, nil
		}

		x1 = x1Next
		x2 = x2Next
	}
}

func main() {
	x1_0 := 6.0
	x2_0 := 1.0

	epsilons := []float64{1e-2, 1e-4, 1e-6, 1e-8}

	for _, eps := range epsilons {
		fmt.Printf("Точность: %v\n", eps)

		root1Newton, root2Newton, itersNewton, errNewton := newtonMethod(x1_0, x2_0, eps)
		if errNewton != nil {
			fmt.Printf("Метод Ньютона: %v\n", errNewton)
		} else {
			fmt.Printf("Метод Ньютона: x1 = %.8f, x2 = %.8f, итераций = %d\n", root1Newton, root2Newton, itersNewton)
		}

		root1Simple, root2Simple, itersSimple, errSimple := simpleIterationMethod(x1_0, x2_0, eps)
		if errSimple != nil {
			fmt.Printf("Метод простой итерации: %v\n\n", errSimple)
		} else {
			fmt.Printf("Метод простой итерации: x1 = %.8f, x2 = %.8f, итераций = %d\n\n", root1Simple, root2Simple, itersSimple)
		}
	}
}