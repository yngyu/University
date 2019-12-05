#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <utility>
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace std;

typedef long long ll;

int main()
{
    int n;
    int i,j;

    double rho = 1.225;
    double T = 293.0;
    double mu = 1.458e-6*pow(T,1.5)/(T+110.4);
    double widthx = 0.5; //単位m
    double u0 = 20.0;
    double Re = rho*u0*widthx/mu;
    double blt = widthx*5.3/sqrt(Re);

    cout << "mu " << mu << endl;
    cout << "Re " << Re << endl;
    
    int nx = 50000;
    int ny = 250;
    double dx = widthx/(double)nx;
    double y_max = 2.0*blt;
    double dy = y_max/(double)ny;
    cout << "y_max " << y_max << endl;
    double dpdx = 0.0; //圧力勾配

    vector<double> y(ny+1);
    for(i = 0;i <= ny;++i){
        y.at(i) = dy*(double)i;
    }

    vector<double> x(nx+1);
    for(i = 0;i <= nx;++i){
        x.at(i) = dx*(double)i;
    }

    vector<vector<double>> u(ny+1,vector<double>(nx+1,0)); //格子
    vector<vector<double>> v(ny+1,vector<double>(nx+1,0)); //格子
    for(i = 0;i <= ny;++i){
        if(i < ny/2){
            u.at(i).at(0) = 0; //はじめの一列初期化
        }else{
            u.at(i).at(0) = u0; //はじめの一列初期化
        }
        v.at(i).at(0) = 0.0; //はじめの一列初期化
    }
    double dudy1;
    double dudy2;
    double dudx;
    double dudx1;
    double dudx2;
    std::chrono::system_clock::time_point  start, end; // 型は auto で可
    start = std::chrono::system_clock::now(); // 計測開始時間
    cout << "for calculate start" << endl;

    for(i = 0;i < nx/2;++i){
        for(j = ny/2;j < ny;++j){
            dudy1 = (u.at(j+1).at(i) - u.at(j-1).at(i))/(2.0*dy);
            dudy2 = (u.at(j+1).at(i) - 2*u.at(j).at(i) + u.at(j-1).at(i))/(dy*dy);
            dudx = (dudy2*mu/rho - v.at(j).at(i)*dudy1 - dpdx/rho)/u.at(j).at(i); //まぁそうか...。
            u.at(j).at(i+1) = u.at(j).at(i) + dx*dudx;
        }
        u.at(ny/2-1).at(i+1) = 0.0; 
        u.at(ny).at(i+1) = u0; //境界条件
        v.at(ny/2-1).at(i+1) = 0;
        for(j = ny/2;j < ny;++j){ //端含まず
            dudx1 = (u.at(j).at(i+1) - u.at(j).at(i))/dx;
            dudx2 = (u.at(j+1).at(i+1) - u.at(j+1).at(i))/dx; //一つ上を取る?y_maxでのvは?不要らしい。
            v.at(j).at(i+1) = v.at(j-1).at(i+1) - (dudx1+dudx2)*dy/2.0; 
        }
        v.at(j).at(i+1) = v.at(j-1).at(i+1); //最終行のみ
        if(i%1000 == 0){
            cout << i << endl;
        }
    }
    for(i = nx/2;i < nx;++i){
       for(j = 1;j < ny;++j){
            dudy1 = (u.at(j+1).at(i) - u.at(j-1).at(i))/(2.0*dy);
            dudy2 = (u.at(j+1).at(i) - 2*u.at(j).at(i) + u.at(j-1).at(i))/(dy*dy);
            if(u.at(j).at(i) != 0){
                dudx = (dudy2*mu/rho - v.at(j).at(i)*dudy1 - dpdx/rho)/u.at(j).at(i); //まぁそうか...。
            }else{
                dudx = 0;
            }
            u.at(j).at(i+1) = u.at(j).at(i) + dx*dudx;
        }

        u.at(0).at(i+1) = 0.0; 
        u.at(ny).at(i+1) = u0; //境界条件
        v.at(0).at(i+1) = 0;

        for(j = 1;j < ny;++j){ //端含まず
            dudx1 = (u.at(j).at(i+1) - u.at(j).at(i))/dx;
            dudx2 = (u.at(j+1).at(i+1) - u.at(j+1).at(i))/dx; //一つ上を取る?y_maxでのvは?不要らしい。
            v.at(j).at(i+1) = v.at(j-1).at(i+1) - (dudx1+dudx2)*dy/2.0; 
        }
        v.at(j).at(i+1) = v.at(j-1).at(i+1); //最終行のみ
        if(i%1000 == 0){
            cout << i << endl;
        } 
    }

    end = std::chrono::system_clock::now();  // 計測終了時間
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count(); //処理に要した時間をミリ秒に変換
    cout << "for caluculate end" << endl;
    cout << elapsed << endl;

    cout << "thickness calculate start" << endl;
    start = std::chrono::system_clock::now();
    /*
    vector<double> boundary_layer_thickness(nx+1);
    vector<double> displacement_thickness(nx+1);
    vector<double> momentum_thickness(nx+1);
    vector<double> energy_thickness(nx+1);
    vector<double> tauw(nx+1);
    vector<double> Cf(nx+1);
    
    for(i = 0;i <= nx;++i){
        for(j = 0;j <= ny;++j){
            if(u.at(j).at(i) >= u0*0.995){
                double temp = j;
                boundary_layer_thickness.at(i) = dy*temp;
                break;
            }
        }
    }
    cout << "boundary end" << endl;
    double temp1;
    double temp2;

    for(i = 0;i <= nx;++i){
        double sum = 0;
        for(j = 0;j < ny;++j){
            temp1 = 1.0-u.at(j).at(i)/u0; //下辺の長さ
            temp2 = 1.0-u.at(j+1).at(i)/u0; //上辺の長さ
            sum += (temp1+temp2)*dy/2.0;
            if(u.at(j+1).at(i) == u0){
                break;
            }
        }
        displacement_thickness.at(i) = sum;
    }
    cout << "displacement end" << endl;

    for(i = 0;i <= nx;++i){
        double sum = 0;
        for(j = 0;j < ny;++j){
            temp1 = u.at(j).at(i)/u0*(1.0-u.at(j).at(i)/u0); //下辺の長さ
            temp2 = u.at(j+1).at(i)/u0*(1.0-u.at(j+1).at(i)/u0); //上辺の長さ
            sum += (temp1+temp2)*dy/2.0;
            if(u.at(j+1).at(i) == u0){
                break;
            }
        }
        momentum_thickness.at(i) = sum;
    }
    cout << "momentum end" << endl;

    for(i = 0;i <= nx;++i){
        double sum = 0;
        for(j = 0;j < ny;++j){
            temp1 = u.at(j).at(i)/u0*(1.0-(u.at(j).at(i)/u0)*(u.at(j).at(i)/u0)); //下辺の長さ
            temp2 = u.at(j+1).at(i)/u0*(1.0-(u.at(j).at(i)/u0)*(u.at(j).at(i)/u0)); //上辺の長さ
            sum += (temp1+temp2)*dy/2.0;
            if(u.at(j+1).at(i) == u0){
                break;
            }
        }
        energy_thickness.at(i) = sum;
    }
    cout << "energy_thickness end" << endl;

    for(i = 0;i <= nx;++i){
        tauw.at(i) = mu*(u.at(1).at(i)-u.at(0).at(i))/dy;
    }
    cout << "tauw end" << endl;
    
    for(i = 0;i <= nx;++i){
        Cf.at(i) = tauw.at(i)/(rho*u0*u0/2.0);
    }
    cout << "Cf end" << endl;

    end = std::chrono::system_clock::now();  // 計測終了時間
    elapsed = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
    cout << elapsed << endl;
    cout << "thickness calculate end" << endl;

    vector<double> boundary_layer_thickness_height(nx+1,0);
    vector<double> displacement_thickness_height(nx+1,0);
    vector<double> momentum_thickness_hieght(nx+1,0);
    vector<double> energy_thickness_height(nx+1,0);
    vector<double> tauw_height(nx+1,0);
    vector<double> Cf_height(nx+1,0);
    for(i = 1;i <= nx;++i){
        boundary_layer_thickness_height.at(i) = boundary_layer_thickness.at(i)/x.at(i)*sqrt(Re/widthx*x.at(i));
        displacement_thickness_height.at(i) = displacement_thickness.at(i)/x.at(i)*sqrt(Re/widthx*x.at(i));
        momentum_thickness_hieght.at(i) = momentum_thickness.at(i)/x.at(i)*sqrt(Re/widthx*x.at(i));
        energy_thickness_height.at(i) = energy_thickness.at(i)/x.at(i)*sqrt(Re/widthx*x.at(i));
        Cf_height.at(i) = Cf.at(i)*sqrt(Re/widthx*x.at(i));
    }
    boundary_layer_thickness_height.at(0) = boundary_layer_thickness_height.at(1);
    displacement_thickness_height.at(0) = displacement_thickness_height.at(1);
    momentum_thickness_hieght.at(0) = momentum_thickness_hieght.at(1);
    energy_thickness_height.at(0) = energy_thickness_height.at(1);
    Cf_height.at(0) = Cf_height.at(1);
    */
    cout << "csv start" << endl;
    
    ofstream ofs("u.csv");
    for(j = 0;j <= ny;++j){
        for(i = 0;i < nx;++i){
            ofs << std::setprecision(50) << u.at(j).at(i) << ",";
        }
        ofs << std::setprecision(50) << u.at(j).at(i) << endl;
        if(j %10 == 0){
            cout << "ofs u:" << j << endl; 
        }
    }

    ofstream ofs2("v.csv");
    for(j = 0;j <= ny;++j){
        for(i = 0;i < nx;++i){
            ofs2 << std::setprecision(50) << v.at(j).at(i) << ",";
        }
        ofs2 << std::setprecision(50) << v.at(j).at(i) << endl;
        if(j %10 == 0){
            cout << "ofs2 v:" << j << endl; 
        }
    }
    /*
    ofstream ofs3("thickness.csv");
    for(i = 0;i < nx;++i){
        ofs3 << std::setprecision(50) << boundary_layer_thickness.at(i) << ",";
        ofs3 << std::setprecision(50) << displacement_thickness.at(i) << ",";
        ofs3 << std::setprecision(50) << momentum_thickness.at(i) << ",";
        ofs3 << std::setprecision(50) << energy_thickness.at(i) << ",";
        ofs3 << std::setprecision(50) << tauw.at(i) << ",";
        ofs3 << std::setprecision(50) << Cf.at(i) << ",";
        ofs3 << std::setprecision(50) << boundary_layer_thickness_height.at(i) << ",";
        ofs3 << std::setprecision(50) << displacement_thickness_height.at(i) << ",";
        ofs3 << std::setprecision(50) << momentum_thickness_hieght.at(i) << ",";
        ofs3 << std::setprecision(50) << energy_thickness_height.at(i) << ",";
        ofs3 << std::setprecision(50) << Cf_height.at(i) << endl;
    }
    ofs3 << std::setprecision(50) << boundary_layer_thickness.at(i) << ",";
    ofs3 << std::setprecision(50) << displacement_thickness.at(i) << ",";
    ofs3 << std::setprecision(50) << momentum_thickness.at(i) << ",";
    ofs3 << std::setprecision(50) << energy_thickness.at(i) << ",";
    ofs3 << std::setprecision(50) << tauw.at(i) << ",";
    ofs3 << std::setprecision(50) << Cf.at(i) << ",";
    ofs3 << std::setprecision(50) << boundary_layer_thickness_height.at(i) << ",";
    ofs3 << std::setprecision(50) << displacement_thickness_height.at(i) << ",";
    ofs3 << std::setprecision(50) << momentum_thickness_hieght.at(i) << ",";
    ofs3 << std::setprecision(50) << energy_thickness_height.at(i) << ",";
    ofs3 << std::setprecision(50) << Cf_height.at(i) << endl;
    cout << "csv end" << endl;
    */
    return 0;
}