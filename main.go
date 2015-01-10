// main
package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"
	"time"
)

const SPARSE_END = 0

type SparseMatrix struct {
	V []float64
	N []int
	C []int
	F []int

	W int
	H int
}

type LUPS struct {
	U  SparseMatrix
	LF []int
	P  []int
	Pt []int
}

func panicIf(err error) {
	if err != nil {
		panic(err)
	}
}

func (m *SparseMatrix) add(r int, c int, v float64) {

	var q = len(m.N)
	m.C = append(m.C, c)
	m.V = append(m.V, v)
	m.N = append(m.N, SPARSE_END)
	if m.F[r] == SPARSE_END {
		m.F[r] = q
	} else {
		m.N[q-1] = q
	}
}

func (m *SparseMatrix) readFromFile(filename string) {

	m.V = append(m.V, 0)
	m.C = append(m.C, 0)
	m.N = append(m.N, SPARSE_END)

	file, err := os.Open(filename)
	panicIf(err)
	defer file.Close()

	reader := bufio.NewReader(file)

	tH, err := reader.ReadString(' ')
	panicIf(err)
	m.H, err = strconv.Atoi(strings.Trim(tH, " "))
	panicIf(err)
	fmt.Println(m.H)

	tW, err := reader.ReadString('\n')
	panicIf(err)
	m.W, err = strconv.Atoi(strings.Trim(tW, "\n"))
	panicIf(err)
	fmt.Println(m.W)

	m.F = make([]int, m.H+1)

	for r := 1; r <= m.H; r += 1 {

		line, err := reader.ReadString('\n')
		panicIf(err)
		pairs := strings.Split(strings.Trim(line, "\n"), "\t")

		for i := 0; i < len(pairs)-1; i += 1 {
			numbers := strings.Split(pairs[i], " ")
			tC, tV := numbers[0], numbers[1]
			c, err := strconv.Atoi(strings.Trim(tC, " "))
			panicIf(err)
			v, err := strconv.ParseFloat(strings.Trim(tV, " "), 64)
			panicIf(err)

			m.add(r, c, v)
		}
	}
}

func (m *SparseMatrix) swapRow(r1, r2 int) {
	m.F[r1], m.F[r2] = m.F[r2], m.F[r1]
}

func max(a, b int) int {
	if a < b {
		return b
	}
	return a
}

func (m *SparseMatrix) swapCol(c1, c2 int) {

	for i := 1; i <= m.H; i += 1 {
		f1, f2 := false, false
		q1, q2 := m.F[i], m.F[i]
		p1, p2 := SPARSE_END, SPARSE_END
		jmax := max(c1, c2)

		q, p := m.F[i], SPARSE_END
		for q != SPARSE_END {
			j := m.C[q]
			if j == c1 {
				q1, p1, f1 = q, p, true
			}
			if !f1 && j < c1 {
				p1, q1 = q, m.N[q]
			}

			if j == c2 {
				q2, p2, f2 = q, p, true
			}
			if !f2 && j < c2 {
				p2, q2 = q, m.N[q]
			}

			if (f1 && f2) || j > jmax {
				break
			}

			q, p = m.N[q], q
		}

		if f1 && f2 {
			m.V[q1], m.V[q2] = m.V[q2], m.V[q1]
		}

		if f1 && !f2 {
			if q1 == q2 || p1 == p2 {
				m.C[q1] = c2
				} else {
				if p1 == SPARSE_END {
					m.F[i] = m.N[q1]
				} else {
					m.N[p1] = m.N[q1]
				}

				if p2 == SPARSE_END {
					m.F[i] = q1
				} else {
					m.N[p2] = q1
				}

				m.N[q1] = q2
				m.C[q1] = c2
			}
		}

		if !f1 && f2 {
			if q1 == q2 || p1 == q2 {
				m.C[q2] = c1
			} else {
				if p2 == SPARSE_END {
					m.F[i] = m.N[q2]
				} else {
					m.N[p2] = m.N[q2]
				}

				if p1 == SPARSE_END {
					m.F[i] = q2
				} else {
					m.N[p1] = q2
				}

				m.N[q2] = q1
				m.C[q2] = c1
			}
		}
	}
}

func LUTriang(m *SparseMatrix, t *LUPS) bool {

	// Матрицы L и U храним в одной _M.U, но _M.U.F будет означать начало U,
	// а _M.LF будет означать начало L.
	if m.H != m.W {
		return false
	}

	H := m.H

	t.U = *m
	t.LF = make([]int, H+1)
	t.P = make([]int, H+H+1)
	t.Pt = make([]int, H+H+1)

	for i := 1; i <= H; i += 1 {
		t.P[i] = i
		t.Pt[H+i] = i
	}

	// Массив "указателей" на последний элемент матрицы L, чтобы можно было
	// в неё переносить элементы из U
	LE := make([]int, H+1)

	Rfill := make([]int, H+1)
	Cfill := make([]int, H+1)

	//Предрассчитаем количество элементов в строках и стообцах
	for k := 1; k <= H; k += 1 {
		for q := t.U.F[k]; q != SPARSE_END; q = t.U.N[q] {
			Rfill[k] += 1
			Cfill[t.U.C[q]] += 1
		}
	}

	Norm := make([]float64, H+1)
	var j int
	var valabs float64

	for k := 1; k <= H; k += 1 {

		//Норма активной подматрицы
		for i := range Norm {
			Norm[i] = 0.0
		}

		for i := k; i <= H; i += 1 {
			for q := t.U.F[i]; q != SPARSE_END; q = t.U.N[q] {
				valabs = math.Abs(t.U.V[q])
				j = t.U.C[q]
				if valabs > Norm[j] {
					Norm[j] = valabs
				}
			}
		}

		optgrowth, opti, optj := H*H, k, k
		optv := 0.0

		for i := k; i <= H; i += 1 {
			for q := t.U.F[i]; q != SPARSE_END; q = t.U.N[q] {
				valabs = math.Abs(t.U.V[q])
				j = t.U.C[q]
				if valabs >= 0.001*Norm[j] {
					growth := (Rfill[i] - 1) * (Cfill[j] - 1)
					if growth < optgrowth || (growth == optgrowth && valabs > optv) {
						optgrowth = growth
						opti = i
						optj = j
						optv = valabs
						if optgrowth == 0 {
							break
						}
					}
				}
			}
			if optgrowth == 0 {
				break
			}
		}

		// Перестановки строк и столбцов
		// при этом перестановки столбцов U не меняют матрицу L
		if opti != k {
			t.U.swapRow(opti, k)
			t.P[opti], t.P[k] = t.P[k], t.P[opti]
			t.LF[opti], t.LF[k] = t.LF[k], t.LF[opti]
			LE[opti], LE[k] = LE[k], LE[opti]
			Rfill[opti], Rfill[k] = Rfill[k], Rfill[opti]
		}
		if optj != k {
			t.U.swapCol(optj, k)
			t.P[H+optj], t.P[H+k] = t.P[H+k], t.P[H+optj]
			Cfill[optj], Cfill[k] = Cfill[k], Cfill[optj]
		}

		var kJ []int
		var kV []float64

		q := t.U.F[k]
		diag := t.U.V[q]

		// Преобразуем строку в U
		// и запомним столбцовые индексы и значения в строке
		for q = t.U.N[q]; q != SPARSE_END; q = t.U.N[q] {
			j := t.U.C[q]
			kJ = append(kJ, j)
			kV = append(kV, t.U.V[q]/diag)

			Cfill[j] -= 1 //Уменьшаем, т.к. убираем элемент в ведущей строке.
		}

		Nz := len(kJ)

		for i := k + 1; i <= H; i += 1 {
			j := t.U.C[t.U.F[i]]
			// Если в строке i нет элемента в веущем столбце, то идем к следующей.
			if j != k {
				continue
			}
			Rfill[i] -= 1 //Уменьшаем, т.к. убираем элемент в ведущем столбце.

			var j1, j2 int
			p := t.U.F[i]
			q := t.U.N[p]

			bdiag := t.U.V[p]
			t.U.V[p] /= diag

			for l := 0; l < Nz; l += 1 {
				j1 = kJ[l]
				for q != SPARSE_END {
					j2 = t.U.C[q]
					if j2 < j1 {
						p, q = q, t.U.N[q]
						continue
					}

					if j2 == j1 {
						t.U.V[q] -= kV[l] * bdiag
					} else {
						// Тут уже нужно добавить элемент
						t.U.C = append(t.U.C, j1)
						p1 := len(t.U.C) - 1
						t.U.N[p] = p1            // предыдущему ссылку на этот
						t.U.N = append(t.U.N, q) // этому ссылку на следующий
						t.U.V = append(t.U.V, -kV[l]*bdiag)
						p, q = q, p1
						// больше вроде ничего делать не нужно,
						// здесь могут попасться только "вторые" элементы в
						// строке и F менять не нужно.

						//Увеличим счетчики элементов
						Rfill[i] += 1
						Cfill[j1] += 1
					}
					break
				}

				if q == SPARSE_END {
					// Добавим элементы в конец
					t.U.C = append(t.U.C, j1)
					p1 := len(t.U.C) - 1
					t.U.N[p] = p1            // предыдущему ссылку на этот
					t.U.N = append(t.U.N, q) // этому ссылку на следующий
					t.U.V = append(t.U.V, -kV[l]*bdiag)
					p, q = q, p1
					// больше вроде ничего делать не нужно,
					// здесь могут попасться только "вторые" элементы в
					// строке и F менять не нужно.

					//Увеличим счетчики элементов
					Rfill[i] += 1
					Cfill[j1] += 1
				}
			}
			// и вот тут элемент в ведущем столбце отходит матрице L
			if t.LF[i] == SPARSE_END {
				t.LF[i] = t.U.F[i]
			} else {
				t.U.N[LE[i]] = t.U.F[i]
			}
			LE[i] = t.U.F[i]
			t.U.F[i] = t.U.N[t.U.F[i]]
			t.U.N[LE[i]] = SPARSE_END
		}
	}

	for i := 1; i <= H; i += 1 {
		t.Pt[t.P[i]] = i
		t.Pt[H+t.P[H+i]] = i
	}

	return true
}

func main() {

	M := new(SparseMatrix)
	M.readFromFile("matrix.txt")

	fmt.Printf("Было ненулевых %v\n", len(M.C)-1)

	T := new(LUPS)

	now := time.Now()

	LUTriang(M, T)

	fmt.Printf("Стало ненулевых %v\n", len(T.U.C)-1)

	elapsed := time.Since(now)
	fmt.Println("Время на разложение ", elapsed.Seconds())

	fmt.Println("ok")
}
