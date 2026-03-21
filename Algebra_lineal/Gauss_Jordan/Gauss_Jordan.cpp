#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

template <typename Tipo> struct Matriz {
  int Renglones, Columnas;
  vector<Tipo> Elementos;
  // Constructor elemento por elemento
  Matriz(int m, int n) : Renglones(m), Columnas(n), Elementos(n * m) {}
  // Constructor por listas
  Matriz(initializer_list<initializer_list<Tipo>> lista) {
    Renglones = lista.size();
    Columnas = lista.begin()->size();
    for (auto &fila : lista) {
      if ((int)fila.size() != Columnas) {
        throw runtime_error("Filas con distinto numero de columnas");
      }
    }

    Elementos.resize(Renglones * Columnas);
    int i = 0;
    for (auto &fila : lista) {
      int j = 0;
      for (auto &xi : fila) {
        (*this)(i, j) = xi;
        j++;
      }
      i++;
    }
  }
  // operador de acceso
  Tipo &operator()(int i, int j) { return Elementos[i * Columnas + j]; }

  Tipo operator()(int i, int j) const { return Elementos[i * Columnas + j]; }
  // operadores algebraicos:suma, resta y multiplicacion
  Matriz operator+(const Matriz &B) const {
    if (Renglones != B.Renglones || Columnas != B.Columnas) {
      throw runtime_error("No puedes sumar matrices de diferentes dimensiones");
    }
    Matriz C(Renglones, Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int j = 0; j < Columnas; j++) {
        C(i, j) = (*this)(i, j) + B(i, j);
      }
    }

    return C;
  }
  Matriz operator-(const Matriz &B) const {
    if (Renglones != B.Renglones || Columnas != B.Columnas) {
      throw runtime_error(
          "No puedes restar matrices de diferentes dimensiones");
    }

    Matriz C(Renglones, Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int j = 0; j < Columnas; j++) {
        C(i, j) = (*this)(i, j) - B(i, j);
      }
    }

    return C;
  }

  Matriz operator*(const Matriz &B) const {
    if (Columnas != B.Renglones) {
      throw runtime_error("No puedes multiplicar estas matrices");
    }
    Matriz C(Renglones, B.Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int j = 0; j < B.Columnas; j++) {
        C(i, j) = 0.0;
        for (int k = 0; k < Columnas; k++) {
          C(i, j) += (*this)(i, k) * B(k, j);
        }
      }
    }
    return C;
  }
};

template <typename Tipo> struct Matriz_Aumentada {
  Matriz<Tipo> &A;
  Matriz<Tipo> B;

  int Renglones;
  int Columnas;
  // Constructor si B es matriz
  Matriz_Aumentada(Matriz<Tipo> &A_, Matriz<Tipo> &B_) : A(A_), B(B_) {
    if (A.Renglones != B.Renglones) {
      throw runtime_error("Numero de filas incompatible");
    }
    Renglones = A.Renglones;
    Columnas = A.Columnas + B.Columnas;
  }
  // Constructor si B es vector
  Matriz_Aumentada(Matriz<Tipo> &A_, const vector<Tipo> &b)
      : A(A_), B(b.size(), 1) {
    if (A.Renglones != (int)b.size()) {
      throw runtime_error("Numeros de filas incompatible");
    }
    for (int i = 0; i < (int)b.size(); i++) {
      B(i, 0) = b[i];
    }
    Renglones = A.Renglones;
    Columnas = A.Columnas + 1;
  }
  // operadores de acceso
  Tipo &operator()(int i, int j) {
    if (j < A.Columnas)
      return A(i, j);
    else
      return B(i, j - A.Columnas);
  }
  const Tipo &operator()(int i, int j) const {
    if (j < A.Columnas)
      return A(i, j);
    else
      return B(i, j - A.Columnas);
  }
};

template <typename Tipo>
vector<Tipo> Gauss_Jordan(Matriz<Tipo> A, vector<Tipo> b) {
  Matriz_Aumentada<Tipo> G(A, b);
  int m = G.Renglones;
  int n = G.Columnas;

  // Pivote
  for (int k = 0; k < m; k++) {
    Tipo pivote = G(k, k);
    // Normalizar fila
    for (int j = 0; j < n; j++) {
      G(k, j) /= pivote;
    }
    // Haciendo ceros
    for (int i = 0; i < m; i++) {
      if (i != k) {
        Tipo factor = G(i, k);
        for (int j = 0; j < n; j++) {
          G(i, j) -= factor * G(k, j);
        }
      }
    }
  }
  vector<Tipo> x(m);
  for (int i = 0; i < m; i++) {
    x[i] = G(i, n - 1);
  }
  return x;
};

int main() {
  Matriz<double> A = {{4, 1, 2}, {2, -2, 2}, {1, -2, 4}};
  vector<double> b = {21, 8, 16};
  vector<double> x = Gauss_Jordan(A, b);
  cout << "La solucion al sistema de ecuaciones es:" << endl;
  for (int i = 0; i < (int)x.size(); i++) {
    cout << "x" << i + 1 << " = " << x[i] << endl;
  }
}
