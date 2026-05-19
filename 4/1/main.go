package main

import (
	"fmt"
	"math"
)

type Point struct {
	X float64
	Y float64
	Z float64
}

// y' = z
func rhsY(x, y, z float64) float64 {
	return z
}

// z' = y'' = 2y + 4x^2*e^(x^2)
func rhsZ(x, y, z float64) float64 {
	return 2.0*y + 4.0*x*x*math.Exp(x*x)
}

func exactY(x float64) float64 {
	s := math.Sqrt(2.0)
	return math.Exp(x*x) + math.Exp(s*x) + math.Exp(-s*x)
}

func buildGrid(a, b, h float64) []float64 {
	n := int(math.Round((b - a) / h))
	xs := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
	}
	return xs
}

func eulerMethod(a, b, h, y0, z0 float64) []Point {
	xs := buildGrid(a, b, h)
	n := len(xs) - 1
	res := make([]Point, n+1)
	res[0] = Point{X: xs[0], Y: y0, Z: z0}

	for k := 0; k < n; k++ {
		xk := res[k].X
		yk := res[k].Y
		zk := res[k].Z

		// y' = z
		// z' = F(x, y, z)
		yNext := yk + h*rhsY(xk, yk, zk)
		zNext := zk + h*rhsZ(xk, yk, zk)

		res[k+1] = Point{
			X: xs[k+1],
			Y: yNext,
			Z: zNext,
		}
	}

	return res
}

func rk4Method(a, b, h, y0, z0 float64) []Point {
	xs := buildGrid(a, b, h)
	n := len(xs) - 1
	res := make([]Point, n+1)
	res[0] = Point{X: xs[0], Y: y0, Z: z0}

	for k := 0; k < n; k++ {
		xk := res[k].X
		yk := res[k].Y
		zk := res[k].Z


		K1 := h * rhsY(xk, yk, zk)
		L1 := h * rhsZ(xk, yk, zk)

		K2 := h * rhsY(xk+h/2.0, yk+K1/2.0, zk+L1/2.0)
		L2 := h * rhsZ(xk+h/2.0, yk+K1/2.0, zk+L1/2.0)

		K3 := h * rhsY(xk+h/2.0, yk+K2/2.0, zk+L2/2.0)
		L3 := h * rhsZ(xk+h/2.0, yk+K2/2.0, zk+L2/2.0)

		K4 := h * rhsY(xk+h, yk+K3, zk+L3)
		L4 := h * rhsZ(xk+h, yk+K3, zk+L3)

		deltaY := (K1 + 2.0*K2 + 2.0*K3 + K4) / 6.0
		deltaZ := (L1 + 2.0*L2 + 2.0*L3 + L4) / 6.0

		res[k+1] = Point{
			X: xs[k+1],
			Y: yk + deltaY,
			Z: zk + deltaZ,
		}
	}

	return res
}

func adams4Method(a, b, h, y0, z0 float64) []Point {
	xs := buildGrid(a, b, h)
	n := len(xs) - 1
	res := make([]Point, n+1)

	startEnd := a + 3.0*h
	if startEnd > b {
		startEnd = b
	}
	start := rk4Method(a, startEnd, h, y0, z0)
	copy(res[:len(start)], start)

	for k := 3; k < n; k++ {
		fk := rhsZ(res[k].X, res[k].Y, res[k].Z)
		fk1 := rhsZ(res[k-1].X, res[k-1].Y, res[k-1].Z)
		fk2 := rhsZ(res[k-2].X, res[k-2].Y, res[k-2].Z)
		fk3 := rhsZ(res[k-3].X, res[k-3].Y, res[k-3].Z)

		yNext := res[k].Y + h*(55.0*res[k].Z-59.0*res[k-1].Z+37.0*res[k-2].Z-9.0*res[k-3].Z)/24.0 // (4.25)

		zNext := res[k].Z + h*(55.0*fk-59.0*fk1+37.0*fk2-9.0*fk3)/24.0

		res[k+1] = Point{
			X: xs[k+1],
			Y: yNext,
			Z: zNext,
		}
	}

	return res
}

func rungeRomberg(fine, coarse []Point, p int) []float64 {
	m := len(coarse)
	rr := make([]float64, m)
	denom := math.Pow(2.0, float64(p)) - 1.0

	for i := 0; i < m; i++ {
		rr[i] = math.Abs(fine[2*i].Y-coarse[i].Y) / denom
	}

	return rr
}

func printTable(title string, values []Point, rr []float64) {
	fmt.Println(title)
	fmt.Printf("%-3s %-8s %-16s %-16s %-16s\n", "k", "x", "y", "|y_exact-y|", "RR")

	maxErr := 0.0
	maxRR := 0.0

	for k, pt := range values {
		err := math.Abs(exactY(pt.X) - pt.Y)
		if err > maxErr {
			maxErr = err
		}

		if k%2 == 0 {
			rrIndex := k / 2
			if rrIndex < len(rr) && rr[rrIndex] > maxRR {
				maxRR = rr[rrIndex]
			}
		}

		if k%2 == 0 {
			rrIndex := k / 2
			if rrIndex < len(rr) {
				fmt.Printf("%-3d %-8.4f %-16.10f %-16.10e %-16.10e\n",
					k, pt.X, pt.Y, err, rr[rrIndex])
			} else {
				fmt.Printf("%-3d %-8.4f %-16.10f %-16.10e %-16s\n",
					k, pt.X, pt.Y, err, "-")
			}
		} else {
			fmt.Printf("%-3d %-8.4f %-16.10f %-16.10e %-16s\n",
				k, pt.X, pt.Y, err, "-")
		}
	}

	fmt.Printf("max exact error = %.10e\n", maxErr)
	fmt.Printf("max RR = %.10e\n\n", maxRR)
}

func main() {
	const (
		a  = 0.0
		b  = 1.0
		h  = 0.1
		y0 = 3.0
		z0 = 0.0
	)

	segments := (b - a) / h
	if math.Abs(segments-math.Round(segments)) > 1e-9 {
		fmt.Println("(b-a)/h must be integer")
		return
	}

	h2 := 2.0 * h
	segments2 := (b - a) / h2
	if math.Abs(segments2-math.Round(segments2)) > 1e-9 {
		fmt.Println("for Runge-Romberg (b-a)/(2h) must also be integer")
		return
	}

	eulerH := eulerMethod(a, b, h, y0, z0)
	euler2H := eulerMethod(a, b, h2, y0, z0)
	rrEuler := rungeRomberg(eulerH, euler2H, 1)

	rk4H := rk4Method(a, b, h, y0, z0)
	rk42H := rk4Method(a, b, h2, y0, z0)
	rrRK4 := rungeRomberg(rk4H, rk42H, 4)

	adamsH := adams4Method(a, b, h, y0, z0)
	adams2H := adams4Method(a, b, h2, y0, z0)
	rrAdams := rungeRomberg(adamsH, adams2H, 4)

	printTable("Euler method", eulerH, rrEuler)
	printTable("Runge-Kutta 4 method", rk4H, rrRK4)
	printTable("Adams 4 method", adamsH, rrAdams)
}