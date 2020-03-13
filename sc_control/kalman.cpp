#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>
#include <vector>
#include <utility>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double I_x = 1.9;
const double I_y = 1.6;
const double I_z = 2.0;

const double diff = 1e-3;

const double omega_s = 17.0*2.0*M_PI/60.0;

void normalize_q(Vector4d& quaternion); //クォータニオン正規化
void true_value_func(vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector); //真値計算&格納
void estimate_value_func(vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector, Matrix<double, 7, 7>& M, Matrix<double, 3, 3>& Q, const Matrix<double, 7, 3>& Gamma, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector); //推定値計算&格納
pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> Runge_kutta(const Vector3d& initial_omega, const Vector4d& initial_q, bool noise_flag); //ルンゲクッタで数値積分
Vector3d calculate_omega_k(const Vector3d& omega, const Vector3d& noise_torque);
double Euler_equation(const double& I_1, const double& I_2, const double& I_3, const double& omega_1, const double& omega_2, const double& noise); //omega計算
Matrix<double, 4,3> calculate_quaternion_Matrix(const Vector4d& quaternion); //クォータニオン計算
void update_M_matrix(Matrix<double, 7, 7>& M, Matrix<double, 3, 3>& Q, const Matrix<double, 7, 3>& Gamma, const pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector); //分散行列アップデート
Matrix<double, 7, 7> Create_A_matrix(const Vector3d& omega, const Vector4d& quaternion);
void store_value(pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value, vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector); //数値格納
Matrix<double, 7, 1> Kalman_Filter(const Vector3d& omega_true, const Vector4d& quaternion_true, const Vector3d& omega_estimate, const Vector4d& quaternion_estimate, Matrix<double, 7, 7>& M, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector); //カルマンフィルター
Vector3d calculate_DCM_vector(const Vector4d& q, const int& flag);
Matrix<double, 3, 7> calculate_H_matrix(const Vector4d& q, const int& flag);
void write_csv(const vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector, const vector<Vector3d, aligned_allocator<Vector3d>>& omega_estimate_store_vector, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_estimate_store_vector, const vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector);
void write_csv_for_omega_and_q(const vector<Vector3d, aligned_allocator<Vector3d>>& omega, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion, bool flag);
void write_csv_for_M(const vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector);

int main()
{
    Vector4d initial_q;
    Vector3d initial_omega;
    initial_q << 1, 0, 0, 0;
    initial_omega << 0.1, omega_s + 0.1, 0;

    int end_time = 50; //end_times[s]秒後に終わる

    int i;

    random_device seed_gen;
    mt19937 mt(seed_gen());
    uniform_real_distribution<> dist(-0.1, 0.1);
    //一様乱数、推定値の初期値用

    Vector4d initial_q_estimate;
    initial_q_estimate << initial_q(0) + dist(mt), initial_q(1) + dist(mt), initial_q(2) + dist(mt), initial_q(3) + dist(mt);
    //initial_q_estimate << initial_q(0) + 0.2, initial_q(1) + 0.2, initial_q(2) + 0.2, initial_q(3) + 0.2;

    normalize_q(initial_q_estimate);
    
    Vector3d initial_omega_estimate;
    initial_omega_estimate << initial_omega(0) + dist(mt), initial_omega(1) + dist(mt), initial_omega(2) + dist(mt);
    //initial_omega_estimate << initial_omega(0) + 0.1, initial_omega(1) + 0.1, initial_omega(2) + 0.1;

    vector<Vector3d, aligned_allocator<Vector3d>> omega_store_vector(1,initial_omega); //真値が入る
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion_store_vector(1,initial_q); //真値が入る

    vector<Vector3d, aligned_allocator<Vector3d>> omega_estimate_store_vector(1,initial_omega_estimate); //推定値が入る
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion_estimate_store_vector(1,initial_q_estimate); //推定値が入る

    //推定分散行列
    Matrix<double, 7, 7> M = MatrixXd::Zero(7, 7);
    for(i = 0;i < 7;++i){
        M(i,i) = 1e-2;
    }

    vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>> M_store_vector(1,M);

    //外乱トルクの分散行列
    Matrix<double, 3, 3> Q = MatrixXd::Zero(3, 3);
    for(i = 0;i < 3;++i){
        Q(i,i) = 1e-4;
    }

    //状態方程式
    Matrix<double, 7, 3> B =  MatrixXd::Zero(7, 3);
    B(4,0) = 1.0/I_x;
    B(5,1) = 1.0/I_y;
    B(6,2) = 1.0/I_z;

    //離散系の行列
    const Matrix<double, 7, 3> Gamma = diff*B;

    for(i = 0;i < end_time;++i){
        //1回あたり、1秒間の計算をする。

        //数値計算
        true_value_func(omega_store_vector, quaternion_store_vector);
        estimate_value_func(omega_estimate_store_vector, quaternion_estimate_store_vector, M, Q, Gamma, M_store_vector);

        Matrix<double, 7, 1> modify;
        modify = Kalman_Filter(omega_store_vector.back(), quaternion_store_vector.back(), omega_estimate_store_vector.back(), quaternion_estimate_store_vector.back(), M, M_store_vector);

        Vector4d modify_quaternion;
        modify_quaternion << modify(0), modify(1), modify(2), modify(3);
        Vector3d modify_omega;
        modify_omega << modify(4), modify(5), modify(6);
        //修正
        quaternion_estimate_store_vector.back() -= modify_quaternion;
        omega_estimate_store_vector.back() -= modify_omega;

        normalize_q(quaternion_estimate_store_vector.back());
    }

    write_csv(omega_store_vector, quaternion_store_vector, omega_estimate_store_vector, quaternion_estimate_store_vector, M_store_vector);

    return 0;
}

void normalize_q(Vector4d& quaternion)
{
    double sum = 0;
    int i;
    for(i = 0;i < 4;++i){
        sum += pow(quaternion(i),2.0);
    }
    sum = sqrt(sum);
    for(i = 0;i < 4;++i){
        quaternion(i) /= sum;
    }
    return;
}

void true_value_func(vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector)
{
    Vector3d initial_omega = omega_store_vector.back();
    Vector4d initial_q = quaternion_store_vector.back();

    //数値格納用
    pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> value;
    //真値なのでnoiseを入れる。
    value = Runge_kutta(initial_omega, initial_q, true);
    //数値格納
    store_value(value, omega_store_vector, quaternion_store_vector);

    return;
}

void estimate_value_func(vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector, Matrix<double, 7, 7>& M, Matrix<double, 3, 3>& Q, const Matrix<double, 7, 3>& Gamma, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector)
{
    Vector3d initial_omega = omega_store_vector.back();
    Vector4d initial_q = quaternion_store_vector.back();

    //数値格納用
    pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> value;
    //推定値、noiseは入れない。
    value = Runge_kutta(initial_omega, initial_q, false);

    //推定分散行列の計算
    update_M_matrix(M, Q, Gamma, value, M_store_vector);

    //数値格納
    store_value(value, omega_store_vector, quaternion_store_vector);
}

pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> Runge_kutta(const Vector3d& initial_omega, const Vector4d& initial_q, bool noise_flag)
{
    int iteration_number = 1.0/diff;

    Vector3d omega;
    omega << initial_omega(0), initial_omega(1), initial_omega(2);
    Vector4d quaternion;
    quaternion << initial_q(0), initial_q(1), initial_q(2), initial_q(3);
    
    int i;

    //random_device seed_gen;
    //mt19937 mt(seed_gen());
    mt19937 mt((int)(omega(0)*1e9)); //MinGWはクソ

    //外乱トルク
    normal_distribution<> dist(0.0, 0.01);
    Vector3d noise_torque;
    noise_torque = Vector3d::Zero();
    
    vector<Vector3d, aligned_allocator<Vector3d>> omega_store_vector;
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion_store_vector;

    for(i = 0;i < iteration_number;++i){
        if(noise_flag){
            noise_torque << dist(mt), dist(mt), dist(mt);
        }

        Vector3d k0 = calculate_omega_k(omega, noise_torque);
        Matrix<double, 4,3> A;
        A = calculate_quaternion_Matrix(quaternion);
        Vector4d l0 = A*omega/2.0*diff;

        Vector3d temp_omega = omega + k0/2.0;
        Vector4d temp_quaternion = quaternion + l0/2.0;

        Vector3d k1 = calculate_omega_k(temp_omega, noise_torque);
        A = calculate_quaternion_Matrix(temp_quaternion);
        Vector4d l1 = A*temp_omega/2.0*diff;

        temp_omega = omega + k1/2.0;
        temp_quaternion = quaternion + l1/2.0;

        Vector3d k2 = calculate_omega_k(temp_omega, noise_torque);
        A = calculate_quaternion_Matrix(temp_quaternion);
        Vector4d l2 = A*temp_omega/2.0*diff;

        temp_omega = omega + k2;
        temp_quaternion = quaternion + l2;

        Vector3d k3 = calculate_omega_k(temp_omega, noise_torque);
        A = calculate_quaternion_Matrix(temp_quaternion);
        Vector4d l3 = A*temp_omega/2.0*diff;

        Vector3d k = (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0;
        Vector4d l = (l0 + 2.0*l1 + 2.0*l2 + l3)/6.0;

        omega += k;
        quaternion += l;
        omega_store_vector.push_back(omega);
        quaternion_store_vector.push_back(quaternion);
        //クォータニオン正規化するか?(しなくてもいいだろ)
    }

    pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>> value = make_pair(omega_store_vector, quaternion_store_vector);
    return value;
}

Vector3d calculate_omega_k(const Vector3d& omega, const Vector3d& noise_torque)
{
    double k_x = Euler_equation(I_y, I_z, I_x, omega(1), omega(2), noise_torque(0));
    double k_y = Euler_equation(I_z, I_x, I_y, omega(2), omega(0), noise_torque(1));
    double k_z = Euler_equation(I_x, I_y, I_z, omega(0), omega(1), noise_torque(2));
    Vector3d k;
    k << k_x, k_y, k_z;
    return k*diff;
}

double Euler_equation(const double& I_1, const double& I_2, const double& I_3, const double& omega_1, const double& omega_2, const double& noise)
{
    double k = (I_1 - I_2)/I_3*omega_1*omega_2 + noise/I_3;
    return k;
}

Matrix<double, 4,3> calculate_quaternion_Matrix(const Vector4d& quaternion)
{
    Matrix<double, 4,3> A;
    A << -quaternion(1), -quaternion(2), -quaternion(3),
        quaternion(0), -quaternion(3), quaternion(2),
        quaternion(3), quaternion(0), -quaternion(1),
        -quaternion(2), quaternion(1), quaternion(0);
    return A;
}

void update_M_matrix(Matrix<double, 7, 7>& M, Matrix<double, 3, 3>& Q, const Matrix<double, 7, 3>& Gamma, const pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector)
{
    vector<Vector3d, aligned_allocator<Vector3d>> omega = value.first;
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion = value.second;
    int n = omega.size(), m = quaternion.size();
    if(n != m){
        cout << "something wrong" << endl;
    }

    int i;

    for(i = 0;i < n;++i){
        Matrix<double, 7, 7> Phi;
        Phi = Create_A_matrix(omega.at(i), quaternion.at(i))*diff;
        Phi += MatrixXd::Identity(7,7);
        M = Phi*M*(Phi.transpose()) + Gamma*Q*(Gamma.transpose());

        M_store_vector.push_back(M); //分散行列保存
    }

    return;
}

Matrix<double, 7, 7> Create_A_matrix(const Vector3d& omega, const Vector4d& quaternion)
{
    Matrix<double, 7, 7> A = MatrixXd::Zero(7, 7);

    A(0,0) = 0; A(0,1) = -0.5*omega(0); A(0,2) = -0.5*omega(1); A(0,3) = -0.5*omega(2); A(0,4) = -0.5*quaternion(1); A(0,5) = -0.5*quaternion(2); A(0,6) = -0.5*quaternion(3);
    A(1,0) = 0.5*omega(0); A(1,1) = 0; A(1,2) = 0.5*omega(2); A(1,3) = -0.5*omega(1); A(1,4) = 0.5*quaternion(0); A(1,5) = -0.5*quaternion(3); A(1,6) = 0.5*quaternion(2);
    A(2,0) = 0.5*omega(1); A(2,1) = -0.5*omega(2); A(2,2) = 0; A(2,3) = 0.5*omega(0); A(2,4) = 0.5*quaternion(3); A(2,5) = 0.5*quaternion(0); A(2,6) = -0.5*quaternion(1);
    A(3,0) = 0.5*omega(2); A(3,1) = 0.5*omega(1); A(3,2) = -0.5*omega(0); A(3,3) = 0; A(3,4) = -0.5*quaternion(2); A(3,5) = 0.5*quaternion(1); A(3,6) = 0.5*quaternion(0);

    A(4,5) = (I_y - I_z)/I_x*omega(2); A(4,6) = (I_y - I_z)/I_x*omega(1);
    A(5,4) = (I_z - I_x)/I_y*omega(2); A(5,6) = (I_z - I_x)/I_y*omega(0);
    A(6,4) = (I_x - I_y)/I_z*omega(1); A(6,5) = (I_x - I_y)/I_z*omega(0);

    return A;
}

void store_value(pair<vector<Vector3d, aligned_allocator<Vector3d>>, vector<Vector4d, aligned_allocator<Vector4d>>>& value, vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector)
{
    vector<Vector3d, aligned_allocator<Vector3d>> omega = value.first;
    vector<Vector4d, aligned_allocator<Vector4d>> quaternion = value.second;
    int n = omega.size();
    int m = quaternion.size();
    if(n != m){
        cout << "something wrong" << endl;
    }

    int i;

    for(i = 0;i < omega.size();++i){
        omega_store_vector.push_back(omega.at(i));
    }
    for(i = 0;i < quaternion.size();++i){
        quaternion_store_vector.push_back(quaternion.at(i));
    }
    return;
}

Matrix<double, 7, 1> Kalman_Filter(const Vector3d& omega_true, const Vector4d& quaternion_true, const Vector3d& omega_estimate, const Vector4d& quaternion_estimate, Matrix<double, 7, 7>& M, vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector)
{
    Vector3d y_true; //観測ベクトル(true)
    Vector3d y_estimate; //観測ベクトル(estimate);

    mt19937 mt((int)(quaternion_true(0)*1e9));
    uniform_int_distribution<> dist_for_flag(0, 2);
    int observe_flag = dist_for_flag(mt); //どの列ベクトルを観測するか

    //観測ノイズ
    Vector3d noise;
    normal_distribution<> dist_for_noise(0.0, 1e-2);
    noise << dist_for_noise(mt), dist_for_noise(mt), dist_for_noise(mt);

    y_true = calculate_DCM_vector(quaternion_true, observe_flag);
    y_estimate = calculate_DCM_vector(quaternion_estimate, observe_flag);
    Vector3d z = y_estimate - y_true;

    z += noise; //ノイズを足す

    Matrix<double, 3, 7> H;
    H = calculate_H_matrix(quaternion_estimate, observe_flag);

    int i;

    Matrix3d R = MatrixXd::Zero(3, 3);
    for(i = 0;i < 3;++i){
        R(i,i) = 1e-4; //観測ノイズの分散
    }
    
    Matrix3d temp;
    temp = H*M*(H.transpose()) + R;

    Matrix<double, 7, 3> K; //カルマンゲイン
    K = M*(H.transpose())*(temp.inverse());

    Matrix<double, 7, 1> modify; //ずれ
    modify = K*z;

    M = (MatrixXd::Identity(7,7) - K*H)*M; //分散更新
    /*
    Matrix3d temp;
    temp = H*M*(H.transpose()) + R;

    Matrix<double, 7, 7> P;
    P = M - M*(H.transpose())*(temp.inverse())*H*M;

    Matrix<double, 7, 3> K;
    K = P*(H.transpose())*(R.inverse());

    Matrix<double, 7, 1> modify; //ずれ
    modify = K*z;

    M = P;
    */
    M_store_vector.back() = M; //保存

    return modify;
}

Vector3d calculate_DCM_vector(const Vector4d& q, const int& flag)
{   
    Vector3d y;
    if(flag == 0){
        y(0) = q(0)*q(0) + q(1)*q(1) - q(2)*q(2) - q(3)*q(3);
        y(1) = 2.0*(q(1)*q(2) + q(0)*q(3));
        y(2) = 2.0*(q(1)*q(3) - q(0)*q(2));
    }else if(flag == 1){
        y(0) = 2.0*(q(1)*q(2) - q(0)*q(3)); 
        y(1) = q(0)*q(0) - q(1)*q(1) + q(2)*q(2) - q(3)*q(3);
        y(2) = 2.0*(q(2)*q(3) + q(0)*q(1));
    }else if(flag == 2){
        y(0) = 2.0*(q(1)*q(3) + q(0)*q(2)); 
        y(1) = 2.0*(q(2)*q(3) - q(0)*q(1));
        y(2) = q(0)*q(0) - q(1)*q(1) - q(2)*q(2) + q(3)*q(3);
    }else{
        cout << "flag something wrong" << endl;
    }
    return y;
}

Matrix<double, 3, 7> calculate_H_matrix(const Vector4d& q, const int& flag)
{
    Matrix<double, 3, 7> H = MatrixXd::Zero(3,7);
    if(flag == 0){
        H(0,0) = q(0); H(0,1) = q(1); H(0,2) = -q(2); H(0,3) = -q(3);
        H(1,0) = q(3); H(1,1) = q(2); H(1,2) = q(1); H(1,3) = q(0);
        H(2,0) = -q(2); H(2,1) = q(3); H(2,2) = -q(0); H(2,3) = q(1);
    }else if(flag == 1){
        H(0,0) = -q(3); H(0,1) = q(2); H(0,2) = q(1); H(0,3) = -q(0);
        H(1,0) = q(0); H(1,1) = -q(1); H(1,2) = q(2); H(1,3) = -q(3);
        H(2,0) = q(1); H(2,1) = q(0); H(2,2) = q(3); H(2,3) = q(2);
    }else if(flag == 2){
        H(0,0) = q(2); H(0,1) = q(3); H(0,2) = q(0); H(0,3) = q(1);
        H(1,0) = -q(1); H(1,1) = -q(0); H(1,2) = q(3); H(1,3) = q(2);
        H(2,0) = q(0); H(2,1) = -q(1); H(2,2) = -q(2); H(2,3) = q(3);
    }else{
        cout << "flag something wrong" << endl;
    }

    H *= 2.0;

    return H;
}

void write_csv(const vector<Vector3d, aligned_allocator<Vector3d>>& omega_store_vector, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_store_vector, const vector<Vector3d, aligned_allocator<Vector3d>>& omega_estimate_store_vector, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion_estimate_store_vector, const vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector)
{
    write_csv_for_omega_and_q(omega_store_vector, quaternion_store_vector, true);
    write_csv_for_omega_and_q(omega_estimate_store_vector, quaternion_estimate_store_vector, false);
    write_csv_for_M(M_store_vector);

    return;
}

void write_csv_for_omega_and_q(const vector<Vector3d, aligned_allocator<Vector3d>>& omega, const vector<Vector4d, aligned_allocator<Vector4d>>& quaternion, bool flag)
{
    int n = omega.size();
    int m = quaternion.size();
    if(n != m){
        cout << "something wrong" << endl;
    }

    string file_name;
    if(flag){
        file_name = "omega_q_true.csv";
    }else{
        file_name = "omega_q_estimate.csv";
    }

    ofstream ofs(file_name);
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

void write_csv_for_M(const vector<Matrix<double, 7, 7>, aligned_allocator<Matrix<double, 7, 7>>>& M_store_vector)
{
    int n = M_store_vector.size();

    ofstream ofs("M_diag.csv");
    int i,j;
    for(i = 0;i < n;++i){
        for(j = 0;j < 6;++j){
            ofs << setprecision(10) << M_store_vector.at(i)(j,j) << ',';
        }
        ofs << setprecision(10) << M_store_vector.at(i)(j,j) << endl;
    }

    return;
}