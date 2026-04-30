package main

import (
	"fmt"
	"math"
)

func makeE(n int) [][]float64 {
	E := make([][]float64, n)
	for i := 0; i < n; i++ {
		E[i] = make([]float64, n)
		E[i][i] = 1.0
	}
	return E
}

func multiply(A, B [][]float64) [][]float64 {
	n := len(A)
	C := make([][]float64, n)
	for i := 0; i < n; i++ {
		C[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				C[i][j] += A[i][k] * B[k][j]
			}
		}
	}
	return C
}

func sign(x float64) float64 {
	if x < 0 {
		return -1.0
	}
	return 1.0
}

func copyMatrix(A [][]float64) [][]float64 {
	n := len(A)
	C := make([][]float64, n)
	for i := 0; i < n; i++ {
		C[i] = make([]float64, len(A[i]))
		copy(C[i], A[i])
	}
	return C
}

func qrDecomposition(A [][]float64) ([][]float64, [][]float64) {
	n := len(A)
	R := copyMatrix(A)
	Q := makeE(n)

	// Выполняем n-1 шагов (по числу столбцов, подлежащих обнулению)
	for k := 0; k < n-1; k++ {
		// Считаем норму ||b||_2 для k-го столбца начиная с элемента k
		sum := 0.0
		for i := k; i < n; i++ {
			sum += R[i][k] * R[i][k]
		}
		normB := math.Sqrt(sum)

		// Строим вектор v по формуле 1.30: v = b + sign(b_1)*||b||*e_1
		v := make([]float64, n)
		for i := k; i < n; i++ {
			v[i] = R[i][k]
		}
		v[k] = v[k] + sign(R[k][k])*normB

		// Считаем v^T * v (скалярное произведение)
		vTv := 0.0
		for i := k; i < n; i++ {
			vTv += v[i] * v[i]
		}

		// Строим матрицу Хаусхолдера H = E - 2*(v*v^T)/(v^T*v)
		H := makeE(n)
		if vTv > 1e-12 { // Защита от деления на ноль, если столбец уже состоит из нулей
			for i := k; i < n; i++ {
				for j := k; j < n; j++ {
					H[i][j] -= 2.0 * v[i] * v[j] / vTv
				}
			}
		}

		// R = H * R
		R = multiply(H, R)
		// Q = Q * H (т.к. H ортогональна и симметрична, транспонирование не нужно)
		Q = multiply(Q, H)
	}
	return Q, R
}

// checkConvergence проверяет критерий остановки для вещественных корней
func checkConvergence(A [][]float64, eps float64) bool {
	n := len(A)
	for m := 0; m < n-1; m++ {
		sum := 0.0
		for l := m + 1; l < n; l++ {
			sum += A[l][m] * A[l][m]
		}
		// Если для какого-то столбца норма поддиагональных элементов больше eps
		// Значит матрица еще не стала треугольной/квазитреугольной
		if math.Sqrt(sum) > eps {
			return false
		}
	}
	return true
}

func main() {
	// Матрица из задания (картинка)
	A_k := [][]float64{
		{5, -5, -6},
		{-1, -8, -5},
		{2, 7, -3},
	}

	eps := 0.01
	n := len(A_k)
	maxIters := 50 // Ограничитель на случай зацикливания из-за комплексных корней

	fmt.Printf("Запуск QR-алгоритма (eps = %v)\n", eps)
	fmt.Println("--------------------------------------------------")

	k := 0
	for {
		// Шаг 1: QR разложение A^{(k)} = Q^{(k)} * R^{(k)}
		Q_k, R_k := qrDecomposition(A_k)

		// Шаг 2: Перемножение в обратном порядке A^{(k+1)} = R^{(k)} * Q^{(k)}
		A_k = multiply(R_k, Q_k)
		k++

		// Проверка критерия останова (упрощенная, для полностью сошедшейся формы)
		// В реальной задаче с комплексными корнями этот критерий может никогда не выполниться
		// для блочной формы 2x2, поэтому мы используем и ограничитель по итерациям.
		if checkConvergence(A_k, eps) || k >= maxIters {
			break
		}
	}

	fmt.Printf("Алгоритм остановлен после %d итераций.\n\n", k)

	fmt.Println("Итоговая матрица A^{(k)}:")
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			fmt.Printf("%8.4f ", A_k[i][j])
		}
		fmt.Println()
	}

	// Извлечение вещественных собственных значений (диагональные элементы)
	fmt.Println("\nСобственные значения (вещественные части лежат на диагонали):")
	for i := 0; i < n; i++ {
		fmt.Printf("λ_%d ≈ %.4f\n", i+1, A_k[i][i])
	}
}