#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace std;

struct All {};
constexpr All all;

template <typename Tipo> struct Vector {
  vector<Tipo> datos;

  Vector(int n) : datos(n) {}

  Vector(initializer_list<Tipo> v) : datos(v) {}

  int size() const { return datos.size(); }
  Tipo &operator[](int i) { return datos[i]; }
  const Tipo &operator[](int i) const { return datos[i]; }

  Vector operator+(const Vector &B) const {
    int n = size(*this);
    int m = size(B);
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
    int n = size(*this);
    int m = size(B);
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

  Vector operator*(Tipo &a) const {
    int n = size(*this);
    Vector C(n);
    for (int i = 0; i < n; i++) {
      C[i] = a * (*this)[i];
    }
    return C;
  };
  Tipo operator*(const Vector &B) const {
    int n = size(*this);
    int m = size(B);
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
};

template <typename Tipo> struct Matriz {
  int Renglones, Columnas;
  Vector<Tipo> Elementos;
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
  // Constructor si B es Vector
  Matriz_Aumentada(Matriz<Tipo> &A_, const Vector<Tipo> &b)
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
template <typename Tipo> struct FVector : public Vector<Tipo> {
  using Vector<Tipo>::Vector;
};
template <typename Tipo> FVector<Tipo> Transpuesta(const Vector<Tipo> &V) {
  FVector R(V.size());
  for (int i = 0; i < V.size(); i++) {
    R[i] = V[i];
  }
  return R;
}
template <typename Tipo> Vector<Tipo> Transpuesta(const FVector<Tipo> &V) {
  FVector R(V.size());
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

template <typename Tipo>
Vector<Tipo> Sustitucion_hacia_atras(const Matriz<Tipo> &U, Vector<Tipo> &y) {
  int n = U.Columnas;
  Vector<Tipo> x(n);
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
Vector<Tipo> Sustitucion_hacia_adelante(const Matriz<Tipo> &L,
                                        Vector<Tipo> &b) {
  int n = L.Columnas;
  Vector<Tipo> x(n);
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
Vector<Tipo> Solver_LU(Matriz<Tipo> &A, Vector<Tipo> &b) {
  int n = A.Renglones;
  Matriz<Tipo> L(n, n);
  Matriz<Tipo> U(n, n);
  tuple<Matriz<Tipo>, Matriz<Tipo>> LUL = LU(A);
  L = get<0>(LUL);
  U = get<1>(LUL);
  Vector<Tipo> y = Sustitucion_hacia_adelante(L, b);
  Vector<Tipo> x = Sustitucion_hacia_atras(U, y);
  return x;
}

template <typename Tipo>
Vector<Tipo> Solver_Gauss_Jordan(Matriz<Tipo> A, Vector<Tipo> b) {
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
  Vector<Tipo> x(m);
  for (int i = 0; i < m; i++) {
    x[i] = G(i, n - 1);
  }
  return x;
};
template <typename Tipo>
Vector<Tipo> Solver_Eliminacion_Gaussiana(Matriz<Tipo> A, Vector<Tipo> b) {
  Matriz_Aumentada<Tipo> G(A, b);
  int n = G.Renglones;
  for (int k = 0; k < n - 1; k++) {
    for (int i = k + 1; i < n; i++) {
      Tipo m = G(i, k) / G(k, k);
      for (int j = k; j <= n; j++) {
        G(i, j) -= m * G(k, j);
      }
    }
  }
  Matriz<Tipo> U = G.A;
  Vector<Tipo> y = G.B;
  return Sustitucion_hacia_atras(U, y);
};
template <typename Tipo> Matriz<Tipo> HouseHolder(const Vector<Tipo> v) {
  int n = v.size();
  FVector<Tipo> w = Transpuesta(v);
  Matriz<Tipo> H(n, n);
  H = Matriz<Tipo>::Identidad(n) - 2.0 * (v * w) / w * v;
  return H;
}

template <typename Tipo>
tuple<Matriz<Tipo>, Matriz<Tipo>> QR(const Matriz<Tipo> A, const Tipo epsilon) {
  int n = A.Renglones;
  Matriz<Tipo> Q(n, n), R(n, n), H(n, n);
  Vector<Tipo> C, e, V;
  R = A;
  Q = Matriz<Tipo>::Identidad(n);
  int i = 0;
  Tipo lambda = Tipo(1);
  Tipo lambda0 = Tipo(0);
  while (lambda - lambda0 > epsilon) {
    lambda0 = lambda;
    C = R(all, i);
    if (C[i] > 0) {
      e = Vector<Tipo>::Zeros(n);
      e[i] = Tipo(1);
    } else if (C[i] < 0) {
      e = Vector<Tipo>::Zeros(n);
      e[i] = Tipo(-1);
    }
    V = C + C.norm() * e;
    H = HouseHolder(V);
    Q = Q * H;
    R = H * R;
    lambda = R(i, i);
    i++;
  }
  return {Q, H};
}

template <typename Tipo>
tuple<vector<Vector<Tipo>>, Vector<Tipo>> Eigen(const Matriz<Tipo> Q,
                                                const Matriz<Tipo> R) {
  int n = Q.Renglones;
  vector<Tipo> lambda(n);
  vector<Vector<Tipo>> v(n);
  for (int i = 0; i < n; i++) {
    lambda[i] = R(i, i);
    v[i] = Q(all, i);
  }
  return {v, lambda};
}
