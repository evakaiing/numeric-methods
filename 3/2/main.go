package main

import (
	"fmt"
)

func main() {
	X := []float64{1.0, 1.9, 2.8, 3.7, 4.6}
	F := []float64{2.4142, 1.0818, 0.50953, 0.11836, -0.24008}
	n := len(X) - 1 // Количество отрезков сплайна (у нас 4)
	x_star := 2.666666667

	// Массивы для коэффициентов сплайна. 
	a := make([]float64, n+1)
	b := make([]float64, n+1)
	c := make([]float64, n+1)
	d := make([]float64, n+1)
	h := make([]float64, n+1)

	for i := 1; i <= n; i++ {
		h[i] = X[i] - X[i-1]
		a[i] = F[i-1]
	}

	alpha := make([]float64, n+1)
	beta := make([]float64, n+1)

	c[1] = 0.0

	for i := 2; i <= n; i++ {
		A := h[i-1]
		B := 2.0 * (h[i-1] + h[i])
		C_coef := 0.0
		if i < n {
			C_coef = h[i]
		}
		
		F_val := 3.0 * ((F[i] - F[i-1])/h[i] - (F[i-1] - F[i-2])/h[i-1])

		alpha[i] = -C_coef / (B + A*alpha[i-1])
		beta[i] = (F_val - A*beta[i-1]) / (B + A*alpha[i-1])
	}

	// обратный
	c[n] = beta[n]
	for i := n - 1; i >= 2; i-- {
		c[i] = alpha[i]*c[i+1] + beta[i]
	}

	for i := 1; i < n; i++ {
		d[i] = (c[i+1] - c[i]) / (3.0 * h[i])
		b[i] = (F[i]-F[i-1])/h[i] - (h[i]*(c[i+1]+2.0*c[i]))/3.0
	}
	d[n] = -c[n] / (3.0 * h[n])
	b[n] = (F[n]-F[n-1])/h[n] - (2.0*h[n]*c[n])/3.0

	fmt.Println("Коэффициенты кубического сплайна")
	fmt.Printf("%-3s %-10s %-10s %-10s %-10s\n", "i", "a_i", "b_i", "c_i", "d_i")
	for i := 1; i <= n; i++ {
		fmt.Printf("%-3d %-10.5f %-10.5f %-10.5f %-10.5f\n", i, a[i], b[i], c[i], d[i])
	}
	fmt.Println()

	fmt.Println("Вычисление в X*:")
	// 3.11
	var interval int = 1
	for i := 1; i <= n; i++ {
		if x_star >= X[i-1] && x_star <= X[i] {
			interval = i
			break
		}
	}
	
	dx_star := x_star - X[interval-1]
	S_star := a[interval] + b[interval]*dx_star + c[interval]*dx_star*dx_star + d[interval]*dx_star*dx_star*dx_star
	
	fmt.Printf("Точка X* = %.9f принадлежит интервалу [%.1f, %.1f]\n", x_star, X[interval-1], X[interval])
	fmt.Printf("Значение сплайна S(X*) = %.6f\n\n", S_star)


	fmt.Println("Точки для графика: ")
	for i := 1; i <= n; i++ {
		step := h[i] / 10.0
		for j := 0; j < 10; j++ { 
			x_val := X[i-1] + float64(j)*step
			dx := x_val - X[i-1]
			y_val := a[i] + b[i]*dx + c[i]*dx*dx + d[i]*dx*dx*dx
			fmt.Printf("%.4f, %.4f\n", x_val, y_val)
		}
	}
	dx_last := h[n]
	y_last := a[n] + b[n]*dx_last + c[n]*dx_last*dx_last + d[n]*dx_last*dx_last*dx_last
	fmt.Printf("%.4f, %.4f\n", X[n], y_last)
}