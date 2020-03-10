#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <utility>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double I_x = 1.9;
const double I_y = 1.6;
const double I_z = 2.0;

pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> Runge_kutta(const double& diff, const int& iteration_number, const Vector3d& initial_omega, const Vector4d& initial_q);
Vector3d calculate_omega_k(const Vector3d& omega, const double& M_c, const double& M_n, const double& diff);
double Euler_equation(const double& I_1, const double& I_2, const double& I_3, const double& omega_1, const double& omega_2, const double& M_c, const double& M_n);
Matrix<double, 4,3> calculate_A(Matrix<double, 4,3>& A, const Vector4d& quaternion);
void write_csv_omega_q(const pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value);

int main()
{
    Vector4d initial_q;
    Vector3d initial_omega;
    const double omega_s = 17.0*2.0*M_PI/60.0;
    initial_q << 1.0, 0.0, 0.0, 0.0;
    initial_omega << 0.1, omega_s + 0.1, 0.0;

    double diff = 1e-2;
    int iteration_number = 5e3;
    pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> value;

    value = Runge_kutta(diff, iteration_number, initial_omega, initial_q);
    write_csv_omega_q(value);

    return 0;
}

pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> Runge_kutta(const double& diff, const int& iteration_number, const Vector3d& initial_omega, const Vector4d& initial_q)
{
    double t_end = diff*iteration_number;
    Vector3d omega;
    omega << initial_omega(0), initial_omega(1), initial_omega(2);
    Vector4d quaternion;
    quaternion << initial_q(0), initial_q(1), initial_q(2), initial_q(3);
    
    int i;

    vector<double> M_c(iteration_number+1,0); //制御トルク、全て0
    vector<double> M_n(iteration_number+1,0); //ノイズトルク、全て0

    vector<Vector3d, aligned_allocator<Vector3d>> omega_store_vector(1,initial_omega);
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion_store_vector(1,initial_q);

    for(i = 0;i < iteration_number;++i){
        Vector3d k0 = calculate_omega_k(omega, M_c.at(i), M_n.at(i), diff);
        Matrix<double, 4,3> A;
        A = calculate_A(A, quaternion);
        Vector4d l0 = A*omega/2.0*diff;

        Vector3d temp_omega = omega + k0/2.0;
        Vector4d temp_quaternion = quaternion + l0/2.0;

        Vector3d k1 = calculate_omega_k(temp_omega, M_c.at(i), M_n.at(i), diff);
        A = calculate_A(A, temp_quaternion);
        Vector4d l1 = A*temp_omega/2.0*diff;

        temp_omega = omega + k1/2.0;
        temp_quaternion = quaternion + l1/2.0;

        Vector3d k2 = calculate_omega_k(temp_omega, M_c.at(i), M_n.at(i), diff);
        A = calculate_A(A, temp_quaternion);
        Vector4d l2 = A*temp_omega/2.0*diff;

        temp_omega = omega + k2;
        temp_quaternion = quaternion + l2;

        Vector3d k3 = calculate_omega_k(temp_omega, M_c.at(i), M_n.at(i), diff);
        A = calculate_A(A, temp_quaternion);
        Vector4d l3 = A*temp_omega/2.0*diff;

        Vector3d k = (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0;
        Vector4d l = (l0 + 2.0*l1 + 2.0*l2 + l3)/6.0;

        omega += k;
        quaternion += l;
        omega_store_vector.push_back(omega);
        quaternion_store_vector.push_back(quaternion);
    }

    pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> value = make_pair(omega_store_vector, quaternion_store_vector);
    return value;
}

Vector3d calculate_omega_k(const Vector3d& omega, const double& M_c, const double& M_n, const double& diff)
{
    double k_x = Euler_equation(I_y, I_z, I_x, omega(1), omega(2), M_c, M_n);
    double k_y = Euler_equation(I_z, I_x, I_y, omega(2), omega(0), M_c, M_n);
    double k_z = Euler_equation(I_x, I_y, I_z, omega(0), omega(1), M_c, M_n);
    Vector3d k;
    k << k_x, k_y, k_z;
    return k*diff;
}

double Euler_equation(const double& I_1, const double& I_2, const double& I_3, const double& omega_1, const double& omega_2, const double& M_c, const double& M_n)
{
    double k = (I_1 - I_2)/I_3*omega_1*omega_2 + (M_c + M_n)/I_3;
    return k;
}

Matrix<double, 4,3> calculate_A(Matrix<double, 4,3>& A, const Vector4d& quaternion)
{
    A << -quaternion(1), -quaternion(2), -quaternion(3),
        quaternion(0), -quaternion(3), quaternion(2),
        quaternion(3), quaternion(0), -quaternion(1),
        -quaternion(2), quaternion(1), quaternion(0);
    return A;
}

void write_csv_omega_q(const pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value)
{
    vector<Vector3d, aligned_allocator<Vector3d>> omega = value.first;
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion = value.second;

    int n = omega.size();
    int m = quaternion.size();
    if(n != m){
        cout << "something wrong" << endl;
    }

    ofstream ofs("omega_q_simulation.csv");
    int i,j;
    for(i = 0;i < 3;++i){
        for(j = 0;j < omega.size()-1;++j){
            ofs << setprecision(10) << omega.at(j)(i) << ",";
        }
        ofs << setprecision(10) << omega.at(j)(i) << endl;
    }

    for(i = 0;i < 4;++i){
        for(j = 0;j < quaternion.size()-1;++j){
            ofs << setprecision(10) << quaternion.at(j)(i) << ",";
        }
        ofs << setprecision(10) << quaternion.at(j)(i) << endl;
    }
    return;
}