package main

import (
	"fmt"
)

func main() {
	X := []float64{-1.0, 0.0, 1.0, 2.0, 3.0}
	Y := []float64{-0.5, 0.0, 0.5, 0.86603, 1.0}
	
	x_target := 2.0
	
	i := 2
	
	xi := X[i]
	xi1 := X[i+1]
	xi2 := X[i+2]
	
	yi := Y[i]
	yi1 := Y[i+1]
	yi2 := Y[i+2]
	
	div1 := (yi1 - yi) / (xi1 - xi)
	
	div2 := (yi2 - yi1) / (xi2 - xi1)
	
	div3 := (div2 - div1) / (xi2 - xi)
	
	y_prime := div1 + div3 * (2.0*x_target - xi - xi1)
	
	y_double_prime := 2.0 * div3
	
	fmt.Printf("y'(%.1f) = %.5f\n", x_target, y_prime)
	fmt.Printf("y''(%.1f) = %.5f\n", x_target, y_double_prime)
}