#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <tuple>
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
    if (lista.size() == 0) {
      throw runtime_error("Lista vacia");
    }
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
  Tipo &operator()(int i, int j) {
    if (i >= Renglones || j >= Columnas) {
      throw out_of_range("Indice fuera de rango");
    }
    return Elementos[i * Columnas + j];
  }

  Tipo operator()(int i, int j) const {
    if (i >= Renglones || j >= Columnas) {
      throw out_of_range("Indice fuera de rango");
    }
    return Elementos[i * Columnas + j];
  }
  static Matriz<Tipo> Zeros(int n, int m) {
    Matriz Z(n, m);
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < m; i++) {
        Z(j, i) = Tipo(0);
      }
    }
    return Z;
  }
  static Matriz<Tipo> Identidad(int n) {
    Matriz I = Zeros(n, n);
    for (int i = 0; i < n; i++) {
      I(i, i) = Tipo(1);
    }
    return I;
  }
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
    Matriz C = Zeros(Renglones, B.Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int k = 0; k < Columnas; k++) {
        for (int j = 0; j < B.Columnas; j++) {
          C(i, j) += (*this)(i, k) * B(k, j);
        }
      }
    }
    return C;
  }
};

template <typename Tipo>
vector<Tipo> Sustitucion_hacia_atras(const Matriz<Tipo> &U, vector<Tipo> &y) {
  int n = U.Columnas;
  vector<Tipo> x(n);
  for (int i = n - 1; i >= 0; i--) {
    x[i] = y[i];
    for (int j = i + 1; j < n; j++) {
      x[i] -= U(i, j) * x[j];
    }
    x[i] /= U(i, i);
  }
  return x;
};
template <typename Tipo>
vector<Tipo> Sustitucion_hacia_adelante(const Matriz<Tipo> &L,
                                        vector<Tipo> &b) {
  int n = L.Columnas;
  vector<Tipo> x(n);
  for (int i = 0; i < n; i++) {
    x[i] = b[i];
    for (int j = 0; j < i; j++) {
      x[i] -= L(i, j) * x[j];
    }
    x[i] /= L(i, i);
  }
  return x;
};
template <typename Tipo> Matriz<Tipo> Eliminacion_Gaussiana(Matriz<Tipo> A) {
  int n = A.Renglones;
  for (int k = 0; k < n - 1; k++) {
    for (int i = k + 1; i < n; i++) {
      Tipo m = A(i, k) / A(k, k);
      // Guardamos los multiplicadores en la parte diagonal inferior de la
      // matriz
      A(i, k) = m;
      for (int j = k + 1; j < n; j++) {
        A(i, j) -= m * A(k, j);
      }
    }
  }
  return A;
};
template <typename Tipo> tuple<Matriz<Tipo>, Matriz<Tipo>> LU(Matriz<Tipo> &A) {
  int n = A.Renglones;
  Matriz<Tipo> L = Matriz<Tipo>::Identidad(n); // Matriz triangular inferior
  Matriz<Tipo> U(n, n);                        // Matriz triangulas superior
  Matriz<Tipo> Gauss = Eliminacion_Gaussiana(A);
  // Separamos las matrices en tringular superior o inferior
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i > j) {
        L(i, j) = Gauss(i, j);
      } else {
        U(i, j) = Gauss(i, j);
      }
    }
  }
  return {L, U};
};

template <typename Tipo>
vector<Tipo> Solver_LU(Matriz<Tipo> &A, vector<Tipo> &b) {
  int n = A.Renglones;
  Matriz<Tipo> L(n, n);
  Matriz<Tipo> U(n, n);
  tuple<Matriz<Tipo>, Matriz<Tipo>> LUL = LU(A);
  L = get<0>(LUL);
  U = get<1>(LUL);
  vector<Tipo> y = Sustitucion_hacia_adelante(L, b);
  vector<Tipo> x = Sustitucion_hacia_atras(U, y);
  return x;
}
int main() {
  Matriz<double> A = {
      {9, -4, -2, 0}, {-4, 17, -6, -3}, {-2, -6, 14, -6}, {0, -3, -6, 11}};
  vector<double> b = {24, -16, 0, 18};
  vector<double> x = Solver_LU(A, b);
  cout << "La solucion al sistema de ecuaciones es:" << endl;
  for (int i = 0; i < (int)x.size(); i++) {
    cout << "x" << i + 1 << " = " << x[i] << endl;
  }
}
