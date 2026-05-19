package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return x / math.Pow(3.0*x+4.0, 3)
}

func integrateRect(X0, Xk, h float64) float64 {
	n := int(math.Round((Xk - X0) / h))
	res := 0.0
	for i := 0; i < n; i++ {
		xMid := X0 + float64(i)*h + h/2.0
		res += f(xMid)
	}
	return res * h
}

func integrateTrap(X0, Xk, h float64) float64 {
	n := int(math.Round((Xk - X0) / h))
	res := 0.0

	for i := 1; i <= n; i++ {
		x_prev := X0 + float64(i-1)*h 
		x_curr := X0 + float64(i)*h   

		res += f(x_prev) + f(x_curr)
	}

	return 0.5 * res * h
}

func integrateSimp(X0, Xk, h float64) float64 {
	n := int(math.Round((Xk - X0) / h))
	if n%2 != 0 {
		fmt.Println("Для метода Симпсона нужно четное число интервалов")
		return math.NaN()
	}

	res := f(X0) + f(Xk)
	for i := 1; i < n; i++ {
		x := X0 + float64(i)*h
		if i%2 == 0 {
			res += 2.0 * f(x)
		} else {
			res += 4.0 * f(x)
		}
	}
	return res * (h / 3.0)
}

func F_exact(x float64) float64 {
	return (-3.0*x - 2.0) / (18.0 * math.Pow(3.0*x+4.0, 2))
}

func main() {
	X0 := -1.0
	Xk := 1.0
	h1 := 0.5
	h2 := 0.25

	rectH1 := integrateRect(X0, Xk, h1)
	trapH1 := integrateTrap(X0, Xk, h1)
	simpH1 := integrateSimp(X0, Xk, h1)

	rectH2 := integrateRect(X0, Xk, h2)
	trapH2 := integrateTrap(X0, Xk, h2)
	simpH2 := integrateSimp(X0, Xk, h2)

	// k = h1 / h2 = 2.0
	k := 2.0
	rungeRect := rectH2 + (rectH2-rectH1)/(math.Pow(k, 2)-1.0)
	rungeTrap := trapH2 + (trapH2-trapH1)/(math.Pow(k, 2)-1.0)
	rungeSimp := simpH2 + (simpH2-simpH1)/(math.Pow(k, 4)-1.0)

	exactVal := F_exact(Xk) - F_exact(X0)

	fmt.Printf("Ззначение: %.5f\n\n", exactVal)

	fmt.Printf("h1 = %.2f\n", h1)
	fmt.Printf("Прямоугольников: %.5f\n", rectH1)
	fmt.Printf("Трапеций: %.5f\n", trapH1)
	fmt.Printf("Симпсона: %.5f\n\n", simpH1)

	fmt.Printf("h2 = %.2f\n", h2)
	fmt.Printf("Прямоугольников: %.5f\n", rectH2)
	fmt.Printf("Трапеций: %.5f\n", trapH2)
	fmt.Printf("Симпсона: %.5f\n\n", simpH2)

	fmt.Printf("Уточнение по Рунге-Ромбергу\n")
	fmt.Printf("Прямоугольников: %.5f (Погрешность: %.5f)\n", rungeRect, math.Abs(exactVal-rungeRect))
	fmt.Printf("Трапеций: %.5f (Погрешность: %.5f)\n", rungeTrap, math.Abs(exactVal-rungeTrap))
	fmt.Printf("Симпсона: %.5f (Погрешность: %.5f)\n", rungeSimp, math.Abs(exactVal-rungeSimp))
}