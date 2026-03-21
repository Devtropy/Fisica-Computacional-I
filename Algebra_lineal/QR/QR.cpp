#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;
// token
struct All {};
constexpr All all;
//=== === === === === Struct Vector === === === ===

template <typename Tipo> struct Vector {
  vector<Tipo> datos;

  Vector(int n) : datos(n) {}

  Vector() = default;

  Vector(initializer_list<Tipo> v) : datos(v) {}

  int size() const { return datos.size(); }
  Tipo &operator[](int i) { return datos[i]; }
  const Tipo &operator[](int i) const { return datos[i]; }

  Vector operator+(const Vector &B) const {
    int n = (*this).size();
    int m = B.size();
    if (n != m) {
      throw runtime_error("No puedes sumar Vectores de diferentes dimensiones");
    }
    Vector C(n);
    for (int i = 0; i < n; i++) {
      C[i] = (*this)[i] + B[i];
    }
    return C;
  };
  Vector operator-(const Vector &B) const {
    int n = (*this).size();
    int m = B.size();
    if (n != m) {
      throw runtime_error(
          "No puedes restar Vectores de diferentes dimensiones");
    }
    Vector C(n);
    for (int i = 0; i < n; i++) {
      C[i] = (*this)[i] - B[i];
    }
    return C;
  };

  Vector operator*(Tipo a) const {
    int n = (*this).size();
    Vector C(n);
    for (int i = 0; i < n; i++) {
      C[i] = a * (*this)[i];
    }
    return C;
  };
  Tipo operator*(const Vector &B) const {
    int n = (*this).size();
    int m = B.size();
    if (n != m) {
      throw runtime_error(" dimensiones diferentes");
    }
    Tipo punto = Tipo(0);
    for (int i = 0; i < n; i++) {
      punto += (*this)[i] * B[i];
    }
    return punto;
  };
  Tipo norm() const { return sqrt((*this) * (*this)); }

  static Vector<Tipo> Zeros(int n) {
    Vector<Tipo> V(n);
    for (int i = 0; i < n; i++) {
      V[i] = Tipo(0);
    }
    return V;
  };

  Vector potencia(const Tipo p) {
    int n = (*this).size();
    Vector C(n);
    for (int i = 0; i < n; i++) {
      C[i] = pow((*this)[i], p);
    }
    return C;
  };
};
template <typename Tipo>
Vector<Tipo> operator*(const Tipo &a, const Vector<Tipo> &v) {
  int n = v.size();
  Vector<Tipo> C(n);

  for (int i = 0; i < n; i++) {
    C[i] = a * v[i];
  }

  return C;
}
//=== === === === === Struct Matriz === === === ===

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
  Vector<Tipo> operator()(All, int j) const {
    if (j >= Columnas)
      throw out_of_range("Columna fuera de rango");

    Vector<Tipo> v(Renglones);
    for (int i = 0; i < Renglones; i++)
      v[i] = (*this)(i, j);

    return v;
  }
  Vector<Tipo> operator()(int i, All) const {
    if (i >= Renglones)
      throw out_of_range("Fila fuera de rango");

    Vector<Tipo> v(Columnas);
    for (int j = 0; j < Columnas; j++)
      v[j] = (*this)(i, j);

    return v;
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
  Matriz operator*(Tipo a) const {

    Matriz C(Renglones, Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int j = 0; j < Columnas; j++) {
        C(i, j) = a * (*this)(i, j);
      }
    }

    return C;
  }
  Matriz operator/(Tipo a) const {

    Matriz C(Renglones, Columnas);

    for (int i = 0; i < Renglones; i++) {
      for (int j = 0; j < Columnas; j++) {
        C(i, j) = (*this)(i, j) / a;
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

  Matriz Inversa() {}
  Tipo Determinante() {}

  Matriz transpuesta() const {
    int n = Renglones;
    int m = Columnas;
    Matriz B(m, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        B(j, i) = (*this)(i, j);
      }
    }
    return B;
  }

  void columnizar(int j, const Vector<Tipo> &v) {
    if (v.datos.size() != Renglones)
      throw runtime_error("Dimension incorrecta");

    for (int i = 0; i < Renglones; i++)
      (*this)(i, j) = v[i];
  }
  void renglonizar(int i, const Vector<Tipo> &v) {
    if (v.datos.size() != Columnas)
      throw runtime_error("Dimension incorrecta");

    if (i >= Renglones)
      throw out_of_range("Renglon fuera de rango");

    for (int j = 0; j < Columnas; j++)
      (*this)(i, j) = v[j];
  }
};
template <typename Tipo>
Matriz<Tipo> operator*(const Tipo a, const Matriz<Tipo> &A) {
  Matriz<Tipo> C(A.Renglones, A.Columnas);

  for (int i = 0; i < A.Renglones; i++)
    for (int j = 0; j < A.Columnas; j++)
      C(i, j) = a * A(i, j);

  return C;
}

//=== === === === === Struct Vector Fila === === === ===

template <typename Tipo> struct FVector : public Vector<Tipo> {
  using Vector<Tipo>::Vector;
};
template <typename Tipo> FVector<Tipo> Transpuesta(const Vector<Tipo> &V) {
  FVector<Tipo> R(V.size());
  for (int i = 0; i < V.size(); i++) {
    R[i] = V[i];
  }
  return R;
}
template <typename Tipo> Vector<Tipo> Transpuesta(const FVector<Tipo> &V) {
  FVector<Tipo> R(V.size());
  for (int i = 0; i < V.size(); i++) {
    R[i] = V[i];
  }
  return R;
}
template <typename Tipo>
Tipo operator*(const FVector<Tipo> &A, const Vector<Tipo> &B) {
  int n = A.size();
  int m = B.size();
  if (n != m) {
    throw runtime_error(" dimensiones diferentes");
  }
  Tipo punto = Tipo(0);
  for (int i = 0; i < n; i++) {
    punto += A[i] * B[i];
  }
  return punto;
};

template <typename Tipo>
Matriz<Tipo> operator*(const Vector<Tipo> &A, const FVector<Tipo> &B) {
  int n = A.size();
  int m = B.size();
  Matriz<Tipo> M(n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      M(i, j) = A[i] * B[j];
    }
  }
  return M;
}

//=== === === === === Algoritmo QR === === === ===

template <typename Tipo> Matriz<Tipo> HouseHolder(Vector<Tipo> v) {
  int n = v.size();
  FVector<Tipo> w = Transpuesta(v);
  return Matriz<Tipo>::Identidad(n) - 2.0 * (v * w) / (w * v);
}
template <typename Tipo>
tuple<Matriz<Tipo>, Matriz<Tipo>> QR(const Matriz<Tipo> &A) {
  int n = A.Renglones;

  Matriz<Tipo> R = A;
  Matriz<Tipo> Q = Matriz<Tipo>::Identidad(n);

  for (int k = 0; k < n - 1; k++) {
    Vector<Tipo> c(n - k);
    for (int i = k; i < n; i++) {
      c[i - k] = R(i, k);
    }
    Tipo beta = (c[0] >= 0 ? 1 : -1) * c.norm();
    Vector<Tipo> e = Vector<Tipo>::Zeros(n - k);
    e[0] = 1;
    Vector<Tipo> v = c + beta * e;
    Vector<Tipo> V = Vector<Tipo>::Zeros(n);
    for (int i = k; i < n; i++)
      V[i] = v[i - k];
    Matriz<Tipo> H = HouseHolder(V);

    R = H * R;
    Q = Q * H;
  }

  return {Q, R};
}

template <typename Tipo>
tuple<Vector<Tipo>, vector<Vector<Tipo>>> QR_Eigen(Matriz<Tipo> A, Tipo epsilon,
                                                   int max_iter = 100000) {
  int n = A.Renglones;
  Matriz<Tipo> V = Matriz<Tipo>::Identidad(n);
  for (int j = 0; j < max_iter; j++) {
    auto [Q, R] = QR(A);

    A = R * Q;
    V = V * Q;

    Tipo error = 0;
    for (int i = 1; i < n; i++) {
      error += abs(A(i, i - 1));
    }
    if (error < epsilon) {
      break;
    }
  }
  Vector<Tipo> lambda(n);
  vector<Vector<Tipo>> v(n);
  for (int i = 0; i < n; i++) {
    lambda[i] = A(i, i);
    v[i] = V(all, i);
  }
  for (int i = 0; i < n - 1; i++) {
    int imin = i;
    for (int j = i + 1; j < n; j++)
      if (lambda[j] < lambda[imin]) {
        imin = j;
      }
    if (imin != i) {
      swap(lambda[i], lambda[imin]);
      swap(v[i], v[imin]);
    }
  }

  return {lambda, v};
}

//=== === === === === Main === === === ===

int main() {
  double c = 3e8; // m/s
  double pi = 355.0 / 113.0;
  double KCH = 3.5652e29; // amu/s^2
  double KCC = 9.5152e29; // amu/s^2
  double mH = 1.0;        // amu
  double mC = 12.0;       // amu
  Matriz<double> k = {{KCH / mH, -KCH / mH, 0.0, 0.0},
                      {-KCH / mC, (KCH + KCC) / mC, -KCC / mC, 0.0},
                      {0.0, -KCC / mC, (KCH + KCC) / mC, -KCH / mC},
                      {0.0, 0.0, -KCH / mH, KCH / mH}};
  int n = k.Renglones;
  auto [w, v] = QR_Eigen(k, 1e-18);
  cout << "Las autovalores son:" << endl;
  for (int i = 0; i < (int)w.size(); i++) {
    cout << "w" << i + 1 << " = " << w[i] << endl;
  }
  cout << "=============================" << endl;
  Vector<double> omega = w.potencia(0.5);
  Vector<double> lambda = omega.potencia(-1) * (2 * pi * c);

  cout << "Las frecuencias son:" << endl;
  for (int i = 0; i < (int)omega.size(); i++) {
    cout << "ω" << i + 1 << " = " << omega[i] << " hz" << endl;
  }
  cout << "=============================" << endl;
  cout << "Las longitudes de onda son:" << endl;
  for (int i = 0; i < (int)lambda.size(); i++) {
    cout << "λ" << i + 1 << " = " << lambda[i] << " m" << endl;
  }
  return 0;
}
