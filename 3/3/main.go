package main

import (
	"fmt"
	"math"
)

func det2(m [2][2]float64) float64 {
	return m[0][0]*m[1][1] - m[0][1]*m[1][0]
}

func det3(m [3][3]float64) float64 {
	return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]) -
		m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
		m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
}

func main() {
	X := []float64{-0.9, 0.0, 0.9, 1.8, 2.7, 3.6}
	Y := []float64{-0.36892, 0.0, 0.36892, 0.85408, 1.7856, 6.3138}
	N := float64(len(X))

	var sx, sx2, sx3, sx4, sy, sxy, sx2y float64
	for i := 0; i < len(X); i++ {
		x := X[i]
		y := Y[i]
		sx += x
		sx2 += x * x
		sx3 += math.Pow(x, 3)
		sx4 += math.Pow(x, 4)
		sy += y
		sxy += x * y
		sx2y += x * x * y
	}

	M1 := [2][2]float64{{N, sx}, {sx, sx2}}
	M1_a0 := [2][2]float64{{sy, sx}, {sxy, sx2}}
	M1_a1 := [2][2]float64{{N, sy}, {sx, sxy}}

	d := det2(M1)
	a0 := det2(M1_a0) / d
	a1 := det2(M1_a1) / d

	fmt.Printf("Коэффициенты: a0 = %.5f, a1 = %.5f\n", a0, a1)
	fmt.Printf("Уравнение: F1(x) = %.5f + %.5f*x\n", a0, a1)

	var Phi1 float64 = 0
	for i := 0; i < len(X); i++ {
		F1_x := a0 + a1*X[i]
		err := F1_x - Y[i]
		Phi1 += err * err
	}
	fmt.Printf("Ф1 = %.5f\n\n", Phi1)


	fmt.Println("Нормальная система:")
	fmt.Printf("%.4f*a0 + %.4f*a1 + %.4f*a2 = %.5f\n", N, sx, sx2, sy)
	fmt.Printf("%.4f*a0 + %.4f*a1 + %.4f*a2 = %.5f\n", sx, sx2, sx3, sxy)
	fmt.Printf("%.4f*a0 + %.4f*a1 + %.4f*a2 = %.5f\n", sx2, sx3, sx4, sx2y)

	M2 := [3][3]float64{
		{N, sx, sx2},
		{sx, sx2, sx3},
		{sx2, sx3, sx4},
	}
	M2_a0 := [3][3]float64{{sy, sx, sx2}, {sxy, sx2, sx3}, {sx2y, sx3, sx4}}
	M2_a1 := [3][3]float64{{N, sy, sx2}, {sx, sxy, sx3}, {sx2, sx2y, sx4}}
	M2_a2 := [3][3]float64{{N, sx, sy}, {sx, sx2, sxy}, {sx2, sx3, sx2y}}

	d3 := det3(M2)
	b0 := det3(M2_a0) / d3
	b1 := det3(M2_a1) / d3
	b2 := det3(M2_a2) / d3

	fmt.Printf("Коэффициенты: a0 = %.5f, a1 = %.5f, a2 = %.5f\n", b0, b1, b2)
	fmt.Printf("Уравнение: F2(x) = %.5f + %.5f*x + %.5f*x^2\n", b0, b1, b2)

	var Phi2 float64 = 0
	for i := 0; i < len(X); i++ {
		F2_x := b0 + b1*X[i] + b2*X[i]*X[i]
		err := F2_x - Y[i]
		Phi2 += err * err
	}
	fmt.Printf("Ф2 = %.5f\n\n", Phi2)
}