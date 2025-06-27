#include <iostream>
#include <fstream>
#include <cmath>
#include </eigen-3.4.0/Eigen/Sparse>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> SpMat; // 疎行列の型
using spmat = Eigen::SparseMatrix <double, Eigen::ColMajor>;
using Solver = Eigen::SparseLU < spmat, Eigen::COLAMDOrdering <int> >;

const int N {500};
const double Delta_x {1.0 / N};

double get_b1(double a, double b) {
    int M = 100; //分割数
    double sum = 0.0;
    double dx = (b - a) / M;
    double x;

    for(int i=0; i < M; i++) {
        x = a + i * Delta_x;
        sum += (x - a) * sin(8.0 * M_PI * x) * dx / Delta_x;
    }

    return sum;
}

double get_b2(double a, double b) {
    int M = 100; //分割数
    double sum = 0.0;
    double dx = (b - a) / M;
    double x;

    for(int i=0; i < M; i++) {
        x = a + i * Delta_x;
        sum += (b - x) * sin(8.0 * M_PI * x) * dx / Delta_x;
    }

    return sum;
}

double integrate(double a, double b) {
    int M = 100; //分割数
    double sum = 0.0;
    double dx = (b - a) / M;
    double x;

    for(int i=0; i < M; i++) {
        x = a + i * Delta_x;
        sum += (b - x) * dx;
    }

    return sum * (162.0 / Delta_x);
}

int main(void) {
    SpMat spmat1(N+1, N+1); // 宣言+サイズの指定
    spmat1.reserve(3 * (N+1)); // メモリの確保

    double x, xp;

    for (int i=1; i < N; i++) {
        x = Delta_x * i;
        xp = Delta_x * (i+1);

        //phi_j
        spmat1.insert(i, i) = -2.0 / Delta_x + integrate(x, xp);

        //phi_j+1
        spmat1.insert(i, i+1) = 1.0 / Delta_x;

        //phi_j-1
        spmat1.insert(i, i-1) = 1.0 / Delta_x + integrate(x, xp);
    }

    spmat1.insert(0, 0) = 1.0;
    spmat1.insert(N, N) = 1.0;

    //右辺ベクトル
    Eigen::VectorXd b(N+1);
    b(0) = 0;
    for(int i=1; i < N+1; i++) {
        b(i) = get_b1(Delta_x * (i-1), Delta_x * i) + get_b2(Delta_x * i, Delta_x * (i+1)); //右辺のベクトル
    }

    Eigen::VectorXd x_vec(N+1); //解を入れる

    /* LU分解 */
    Solver sol;
    sol.analyzePattern( spmat1 );
    sol.factorize( spmat1 );

    //解く
    x_vec = sol.solve( b );

    //ファイル出力
    std::ofstream ofs("output.dat");
    ofs << 0 << " " << x_vec(0) << " " << 0 << std::endl;
    for(int i=1; i < N-1; i++) {
        ofs << Delta_x * i << " " << x_vec(i) << " " << -sin(8.0 * M_PI * Delta_x * i) / (64.0 * M_PI * M_PI - 81.0) << std::endl;
    }
    ofs << 1 << " " << x_vec(N-1) + 1 << " " << 0 << std::endl;
    ofs.close();

    return 0;
}