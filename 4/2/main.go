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

type Row struct {
	K         int
	X         float64
	Y         float64
	Exact     float64
	ExactErr  float64
	RR        *float64
}

func exactY(x float64) float64 {
	return 1.0 + 1.0/x
}

func exactZ(x float64) float64 {
	return -1.0 / (x * x)
}

func buildGrid(a, b, h float64) []float64 {
	n := int(math.Round((b - a) / h))
	xs := make([]float64, n+1)
	for i := 0; i <= n; i++ {
		xs[i] = a + float64(i)*h
	}
	return xs
}

// y' = z
func rhsY(x, y, z float64) float64 {
	return z
}

// z' = y'' = 2y / (x^2 (x+1))
func rhsZ(x, y, z float64) float64 {
	return 2.0 * y / (x * x * (x + 1.0))
}

func rk4System(a, b, h, y0, z0 float64) []Point {
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

// 2*y(2) - 4*y'(2) - 4 = 0
func shootingResidual(h, eta float64) (float64, []Point) {
	sol := rk4System(1.0, 2.0, h, eta, -1.0)
	last := sol[len(sol)-1]
	phi := 2.0*last.Y - 4.0*last.Z - 4.0
	return phi, sol
}

// y' = z
// z' = -p(x)z - q(x)y + f(x)
func shootingMethod(h, eta0, eta1, tol float64, maxIter int) ([]Point, float64, []float64, []float64) {
	var historyEta []float64
	var historyPhi []float64

	phi0, _ := shootingResidual(h, eta0)
	phi1, sol1 := shootingResidual(h, eta1)

	historyEta = append(historyEta, eta0, eta1)
	historyPhi = append(historyPhi, phi0, phi1)

	if math.Abs(phi0) < tol {
		sol0 := rk4System(1.0, 2.0, h, eta0, -1.0)
		return sol0, eta0, historyEta, historyPhi
	}
	if math.Abs(phi1) < tol {
		return sol1, eta1, historyEta, historyPhi
	}

	for iter := 0; iter < maxIter; iter++ {
		denom := phi1 - phi0
		if math.Abs(denom) < 1e-15 {
			break
		}

		eta2 := eta1 - (eta1-eta0)*phi1/denom

		phi2, sol2 := shootingResidual(h, eta2)

		historyEta = append(historyEta, eta2)
		historyPhi = append(historyPhi, phi2)

		if math.Abs(phi2) < tol {
			return sol2, eta2, historyEta, historyPhi
		}

		eta0, phi0 = eta1, phi1
		eta1, phi1 = eta2, phi2
		sol1 = sol2
	}

	return sol1, eta1, historyEta, historyPhi
}


func finiteDifferenceMethod(h float64) []Point {
	a, b := 1.0, 2.0
	xs := buildGrid(a, b, h)
	n := len(xs) - 1

	low := make([]float64, n+1) // a_i
	diag := make([]float64, n+1) // b_i
	up := make([]float64, n+1) // c_i
	rhs := make([]float64, n+1) // d_i

	diag[0] = -1.0
	up[0] = 1.0
	rhs[0] = -h

	for k := 1; k <= n-1; k++ {
		xk := xs[k]
		q := -2.0 / (xk * xk * (xk + 1.0))

		low[k] = 1.0
		diag[k] = -2.0 + h*h*q
		up[k] = 1.0
		rhs[k] = 0.0
	}


	// y'(2) = (yN - yN-1)/h

	// тк 2y(2) - 4y'(2) = 4
	// 2yN - 4*(yN - yN-1)/h = 4
	low[n] = 4.0 / h
	diag[n] = 2.0 - 4.0/h
	up[n] = 0.0
	rhs[n] = 4.0

	ys := solveTridiagonal(low, diag, up, rhs)

	res := make([]Point, n+1)
	for i := 0; i <= n; i++ {
		z := math.NaN()
		if i == 0 {
			z = (ys[1] - ys[0]) / h
		} else if i == n {
			z = (ys[n] - ys[n-1]) / h
		} else {
			z = (ys[i+1] - ys[i-1]) / (2.0 * h)
		}
		res[i] = Point{X: xs[i], Y: ys[i], Z: z}
	}

	return res
}


func solveTridiagonal(a, b, c, d []float64) []float64 {
	n := len(d)
	cp := make([]float64, n)
	dp := make([]float64, n)

	cp[0] = c[0] / b[0]
	dp[0] = d[0] / b[0]

	for i := 1; i < n; i++ {
		den := b[i] - a[i]*cp[i-1]
		if i < n-1 {
			cp[i] = c[i] / den
		}
		dp[i] = (d[i] - a[i]*dp[i-1]) / den
	}

	x := make([]float64, n)
	x[n-1] = dp[n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = dp[i] - cp[i]*x[i+1]
	}
	return x
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

func buildRows(sol []Point, rr []float64) []Row {
	rows := make([]Row, len(sol))
	for k, pt := range sol {
		err := math.Abs(exactY(pt.X) - pt.Y)

		var rrPtr *float64
		if k%2 == 0 {
			idx := k / 2
			if idx < len(rr) {
				val := rr[idx]
				rrPtr = &val
			}
		}

		rows[k] = Row{
			K:        k,
			X:        pt.X,
			Y:        pt.Y,
			Exact:    exactY(pt.X),
			ExactErr: err,
			RR:       rrPtr,
		}
	}
	return rows
}

func printRows(title string, rows []Row) {
	fmt.Println(title)
	fmt.Printf("%-3s %-8s %-16s %-16s %-16s %-16s\n", "k", "x", "y", "y_exact", "|err|", "RR")

	maxErr := 0.0
	maxRR := 0.0

	for _, r := range rows {
		if r.ExactErr > maxErr {
			maxErr = r.ExactErr
		}
		if r.RR != nil && *r.RR > maxRR {
			maxRR = *r.RR
		}

		rrStr := "-"
		if r.RR != nil {
			rrStr = fmt.Sprintf("%.10e", *r.RR)
		}

		fmt.Printf("%-3d %-8.4f %-16.10f %-16.10f %-16.10e %-16s\n",
			r.K, r.X, r.Y, r.Exact, r.ExactErr, rrStr)
	}

	fmt.Printf("max exact error = %.10e\n", maxErr)
	fmt.Printf("max RR = %.10e\n\n", maxRR)
}

func printShootingIterations(etas, phis []float64) {
	fmt.Println("Shooting iterations")
	fmt.Printf("%-3s %-16s %-16s\n", "j", "eta", "Phi(eta)")
	for i := 0; i < len(etas); i++ {
		fmt.Printf("%-3d %-16.10f %-16.10e\n", i, etas[i], phis[i])
	}
	fmt.Println()
}

func main() {
	const (
		h   = 0.1
		tol = 1e-10
	)

	shootH, eta, etas, phis := shootingMethod(h, 1.5, 2.5, tol, 50)
	shoot2H, _, _, _ := shootingMethod(2.0*h, 1.5, 2.5, tol, 50)

	rrShoot := rungeRomberg(shootH, shoot2H, 4)
	rowsShoot := buildRows(shootH, rrShoot)

	fmt.Printf("Exact y(1) = %.10f\n", exactY(1.0))
	fmt.Printf("Found eta  = %.10f\n\n", eta)
	printShootingIterations(etas, phis)
	printRows("Shooting method", rowsShoot)


	fdH := finiteDifferenceMethod(h)
	fd2H := finiteDifferenceMethod(2.0 * h)

	rrFD := rungeRomberg(fdH, fd2H, 1)
	rowsFD := buildRows(fdH, rrFD)

	printRows("Finite difference method", rowsFD)
}