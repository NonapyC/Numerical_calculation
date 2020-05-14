// ---------------------------------------------------------------------------------------------------
// NONAPY VERSION 2.7.0
// Latest Edit: 18/12/23.
//
// Copyright (C) 2018 Nonaka Soshi †
// † Master of Science, Department of Physics, Osaka University, Japan.
//
// This code is library for algebraic computation.  ( Library for linear algebra is NSLin.h )
//
// Pleare e_mail me for any inquiry : sou.235711_at_gmail.com ( exchange: _at_ -> @ )
// And the first time I corded this is 2017/12/19.
//---------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>

namespace LibNS2{
#include "NSLin.h"
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// 　　                                     PART1. 代数

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ------------------------------------------------------------------------------------------------
//    sqare of x ; sqr(x) = x^2　
// ------------------------------------------------------------------------------------------------
class sqr {
public:
    template<typename Tin>
    inline Tin operator()(const Tin x) const{
        return x * x;
    }
}sqr;

// ------------------------------------------------------------------------------------------------
//    cubic of x ; cubed(x) = x^3 　
// ------------------------------------------------------------------------------------------------
class cube {
public:
    template<typename Tin>
    inline Tin operator()(const Tin x) const{
        return x * x * x;
    }
}cube;

// ---------------------------------------------------------------------------------------------------
//    Gauss-Legendre integration:
//
//    Return the integral of the function between [a,b]
//    by 64-point Gauss-Legendre method.
//    This function is evaluated exactly 64 times at interior points in the range of integration.
//
//    qgaus(func, a, b, args...=0):
//         --->parameter: func : a function which we find a integral of.
//                        a,b : the area [a,b] where we integrate in.
//                        args : (optional) extra arguments passed to the objective function.
//    example -> cd ../../program/cpp/test_integ.cpp
// ---------------------------------------------------------------------------------------------------

template<typename TFunc , typename Tin , typename Tout , typename... Args>
inline Tout qgaus(TFunc func, const Tin a, const Tout b, const Args... args){
    static double ax[] = {0.0,0.024350293,0.072993122,0.121462819,0.16964442,0.217423644,0.264687162,0.311322872,0.357220158,0.402270158,0.446366017,0.489403146,0.531279464,0.571895646,0.611155355,0.648965471,0.685236313,0.71988185,0.752819907,0.783972359,0.813265315,0.840629296,0.865999398,0.889315446,0.910522137,0.929569172,0.946411375,0.9610088,0.973326828,0.983336254,0.991013371,0.996340117,0.999305042};
    static double w[] = {0.0,0.048690957,0.048575467,0.048344762,0.047999389,0.047540166,0.046968183,0.046284797,0.045491628,0.044590558,0.043583725,0.042473515,0.041262563,0.039953741,0.038550153,0.037055129,0.035472213,0.033805162,0.032057928,0.030234657,0.028339673,0.02637747,0.024352703,0.022270174,0.020134823,0.017951716,0.01572603,0.013463048,0.011168139,0.00884676,0.006504458,0.004147033,0.001783281};

    const Tin xm = 0.5 * ( b + a );
    const Tin xr = 0.5 * ( b - a );
    Tin s = 0.0;
    for (int h = 1 ; h <= 32 ; h++ ) {
        Tin dx = xr * ax[h];
        s += w[h] * ( func( xm + dx , args... ) + func( xm - dx , args... ) );
    }
    return s *= xr; /*Scale the answer to the range of integration.*/
}

// ---------------------------------------------------------------------------------------------------
// 　Solve the nonlinear equation by newton method.
//   newton(f, init, ans, args...) : Find a root of a function.
//      ---> parameter ; f : a 1-variable function which we find a root of.
//                       init : initial guess.
//                       ans : a reference where we put the root in.
//                       args : extra arguments passed to the objective function.
//  --example--
//  cd ../../program/cpp/test_fnewton.cpp
// ---------------------------------------------------------------------------------------------------
template<typename TFunc, typename Tin, typename... Args>
inline Tin newton(TFunc func, const Tin init, Tin &ans, const Args... args){
    Tin x_old, x_new;
    const double eps = 1.0e-8;
    const double dx = 1.0e-5;
    x_new = init;
    do {
        x_old = x_new;
        x_new = x_old - func( x_old , args... ) / ( ( func( x_old + dx , args... ) - func( x_old , args... ) ) / dx );
    } while( fabs( x_new - x_old ) > eps || fabs( func( x_new , args... ) ) > eps );
    ans = x_new;
    return 0;
}

// ---------------------------------------------------------------------------------------------------
// 　Solve the nonlinear equation by bisection method
//   R_nibun(f, a_init, b_init, ans, args...) : Find a root of a function.
//      ---> parameter ; f : a 1-variable function to find a root of.
//                       a_init, b_init : initial guess, a < b.
//                       ans : a reference where we put the root in.
//                       args : extra arguments passed to the objective function.
//  --example--
//  cd ../../program/cpp/test_fnewton.cpp
// ---------------------------------------------------------------------------------------------------
template<typename TFunc, typename Tin, typename... Args>
inline Tin R_nibun(TFunc func, Tin a_init, Tin b_init, double &ans, const Args... args){
    const double eps = 1.0e-8;
    if( func( a_init , args... ) * func( b_init , args... ) > 0.0 ) { cout << "Error!! Initial value must be across the solution." << endl; return 0;}
    else{
        while( fabs( a_init - b_init ) > eps ){
            double c = ( a_init + b_init ) / 2.0;
            if( func( a_init , args... ) * func( c , args... ) < 0.0 ) b_init = c;
            else a_init = c;
        }
        ans = a_init;
    }
    return 0;
}

// --------------------------------------------------------------------------------------------------------
// 　2変数連立方程式をNewton法で解く
//
//  fnewton(f, g, x_initial, vy_initial, &ans1, &ans2, args_x, args_y, epsilon) : Find a root of a vector function.
//      ---> parameter ; f,g : 2-variable functions to find a root of. //x(y)_initial : initial guess.
//                       &ans1(2) : reference where we put the root in. //epsilon : (optional) Tolerance for termination.
//                       args_x(y) : a vector, Extra arguments passed to the objective function.
//  --example--
//  cd ../../program/cpp/test_fnewton.cpp
// --------------------------------------------------------------------------------------------------------
template<typename TFunc, typename Tin, typename Tout, typename... Args>
inline Tout xpartial(TFunc func, Tin x, Tout y, const vector<Args...>& arg1, const double epsilon = 1.0e-3){
    return ( func( x + epsilon , y , arg1 ) - func( x , y , arg1 ) ) / epsilon;
}
template<typename TFunc, typename Tin, typename Tout, typename... Args>
inline Tout ypartial(TFunc func, Tin x, Tout y, const vector<Args...>& arg2, const double epsilon = 1.0e-3){
    return ( func( x , y + epsilon , arg2 ) - func( x , y , arg2 ) ) / epsilon;
}
template<typename TFunc_f, typename TFunc_g , typename Tin, typename Tout, typename... Args1, typename... Args2>
inline Tout fnewton(TFunc_f f, TFunc_g g, const Tin x_initial, const Tout y_initial, Tout &ans1, Tout &ans2, const vector<Args1...>& args_x, const vector<Args2...>& args_y, const double epsilon = 1.0e-10){
    Tin x_old, y_old;
    Tin x_new = x_initial;
    Tin y_new = y_initial;
    do {
        x_old = x_new;
        y_old = y_new;
        x_new = x_old + ( f( x_old , y_old, args_x ) * ypartial( g , x_old , y_old , args_y ) - g( x_old ,  y_old , args_y) * ypartial( f , x_old , y_old , args_x ) ) / ( ypartial( f , x_old , y_old , args_x ) * xpartial( g , x_old , y_old , args_y ) - xpartial( f , x_old , y_old , args_x ) * ypartial( g , x_old , y_old , args_y ) );
        y_new = y_old + ( g( x_old , y_old , args_y ) * xpartial( f , x_old , y_old , args_x ) - f( x_old , y_old , args_x ) * xpartial( g , x_old , y_old , args_y ) ) / ( ypartial( f , x_old , y_old , args_x ) * xpartial( g , x_old , y_old , args_y ) - xpartial( f , x_old , y_old , args_x ) * ypartial( g , x_old , y_old , args_y ) );
    } while( hypot( y_new - y_old , x_new - x_old ) > epsilon );

    ans1 = x_new;
    ans2 = y_new;
    return 0;
}

// --------------------------------------------------------------------------------------------------------
// 　N変数 N元 非線形連立方程式をNewton法で解く
//  f_i(x_1, x_2, ..., x_N) = 0 の N変数 N個の非線形連立方程式の解を求める
//  サンプルは、cd ../../program/cpp/test_Nfnewton.cpp
// --------------------------------------------------------------------------------------------------------
template<typename TFunc, typename Tout> 
inline Tout jac(TFunc func, const vector<Tout> n, const int i, const int j){
    const double dn = 1.0e-8;
    vector<double> n_adv(n);
    n_adv[j] += dn;
    return ( func( n_adv , i ) - func( n , i ) ) / dn;
}

template<typename TFunc, typename Tout, typename Tin>
inline Tout Nfnewton(TFunc func, const vector<Tin> init, vector<Tout> &ans, const int Ndim){
    vector<Tin> x_new(init);
    const double eps = 1.0e-8;
    double count = 0.0;
    do {
        vector<vector<double>> J( Ndim, vector<double>( Ndim, 0.0 ) );
        vector<double> b( Ndim , 0.0 );
        vector<double> dn( Ndim , 0.0 );

        count = 0.0;
        vector<Tin> x_old(x_new);
        for (int i = 0; i < Ndim; i++) {
            b[i] = func(x_new, i);
            for (int j = 0; j < Ndim; j++) {
                J[i][j] = jac( func , x_new , i , j );
            }
        }
        LinNS2::LUdecomp( J , Ndim );
        LinNS2::LUSolve( J , b , dn , Ndim );
        for (int i = 0; i < Ndim; i++) {
            x_new[i] = x_old[i] - dn[i];
            count += fabs(dn[i]);
        }
    }while( count > eps );
    for (int i = 0; i < Ndim; i++) {
        ans[i] = x_new[i];
    }
    return 0;
}

// ---------------------------------------------------------------------------------------------------
// 　微分方程式を （4段4次の）Runge-Kutta法 で解く
//
//   dy(t)/dt = f1(y(t),t) の形の常微分方程式を解くのが RK4_1.
//   dy(t)/dt = f1(y(t), z(t),t), dz(t)/dt = f2(y(t),z(t),t). の形の連立常微分方程式を解くのが RK4_2.
//
//   RK4_1(f1, t, y_init, ref_time, ref_y, t_init, eps)
//      ---> parameter ; f1 : dy(t)/dt = f1(y(t),t).  //t : a maximam time at which we stop update.
//                       y_init : initial condition y(0).　//ref_time, ref_y : time array and root array y(t).
//                       t_init : (optiomal) initial time t_0 = 0.0.　//eps : (optional) Tolerance for termination of time.
//
//   RK4_2(f1,f2, t, y_init, z_init, ref_time, ref_y, ref_z, t_init, eps)
//        ---> parameter ; f1 : dy(t)/dt = f1(y(t),z(t),t). //f2 : dz(t)/dt = f2(y(t),z(t),t).
//                       t : a maximam time at which we stop update. //y_init ,z_init: initial condition y(0) and z(0).
//                       ref_time, ref_y, ref_z: time array and root array y(t). //t_init : (optiomal) initial time t_0 = 0.0.
//                       eps : (optional) Tolerance for termination of time.
//
//   --example--
//   double f1(double y, double z, double t){return z;} double f2(double y, double z, double t){return - 1.0 * y - 0.2 * z ;}
//   int main(){ double ref1[10000], ref2[10000] , ref3[10000]; RK4_2( f1, f2, 100. , 0.3 , 0.1 , ref1, ref2, ref3 );}
// ---------------------------------------------------------------------------------------------------
//一階微分のみの微分方程式
template <typename TFunc, typename Tout, typename Tin>
inline void RK4_1(TFunc f1, Tin t, Tout y_init, Tout *ref_time , Tout *ref_y , Tin t_init = 0.0 , const double eps = 1.0e-2){
    Tin k1, k2, k3, k4;
    for (int i = 0; i < int( t / eps ) ; i++) {
        ref_time[i] = t_init;
        ref_y[i] = y_init;
        k1 = eps * f1( y_init , t_init );
        k2 = eps * f1( y_init + k1 / 2. , t_init + eps / 2. );
        k3 = eps * f1( y_init + k2 / 2. , t_init + eps / 2. );
        k4 = eps * f1( y_init + k3 , t_init + eps );
        y_init = y_init + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6. ;
        t_init = t_init + eps;
    }
}
//2階微分を含んだ微分方程式
template <typename TFunc, typename Tout, typename Tin>
inline void RK4_2(TFunc f1, TFunc f2, Tin t, Tin y_init, Tout z_init, Tout *ref_time , Tout *ref_y, Tout *ref_z , Tin t_init = 0.0 , const double eps = 1.0e-2){
    Tin k1, k2, k3, k4, p1, p2, p3, p4;
    for (int i = 0; i < int( t / eps ) ; i++) {
        ref_time[i] = t_init;
        ref_y[i] = y_init; 
        ref_z[i] = z_init;
        k1 = eps * f1( y_init , z_init , t_init );
        p1 = eps * f2( y_init , z_init , t_init );

        k2 = eps * f1( y_init + k1 / 2. , z_init + p1 / 2. , t_init + eps / 2. );
        p2 = eps * f2( y_init + k1 / 2. , z_init + p1 / 2. , t_init + eps / 2. );

        k3 = eps * f1( y_init + k2 / 2. , z_init + p2 / 2. , t_init + eps / 2. );
        p3 = eps * f2( y_init + k2 / 2. , z_init + p2 / 2., t_init + eps / 2. );

        k4 = eps * f1( y_init + k3 , z_init + p3 , t_init + eps );
        p4 = eps * f2( y_init + k3 , z_init + p3 , t_init + eps );

        y_init = y_init + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6. ;
        z_init = z_init + ( p1 + 2.0 * p2 + 2.0 * p3 + p4 ) / 6. ;
        t_init = t_init + eps;
    }
}

//---------------------------------------------------------------------------------------------
//    離散フーリエ変換：
//          x空間の関数f(x)を波数k(=2π/λ)空間の関数F(k)に変換する。
//
//    ftfメソッド：関数をdxの幅のnumber_of_x_sampling個の離散データに離散化した後に離散フーリエ変換
//    ftfataメソッド：幅dxをもつnumber_of_x_sampling個の離散データを離散フーリエ変換
//    サンプルは、cd ../../program/cpp/test_FT.cpp
//
//---------------------------------------------------------------------------------------------
template <typename TFunc>//関数をdxの幅で離散化する
void f_descrete(TFunc func, vector<double>& vec, const double dx, const int number_of_x_sampling){
    for (int i = 0; i < number_of_x_sampling; i++) {
        vec[i] = func( i * dx );
    }
}

class DFT{//離散フーリエ変換用のクラス
private:
    vector<double> vec_x;
public:
    vector<double> ftReal;
    vector<double> ftImag;
    vector<double> ftAbs;

    template <typename TFunc>
    inline double ftf(TFunc func, const double dx, const int number_of_x_sampling){
        for (int i = 0; i < number_of_x_sampling; i++) {//コンストラクタ
            ftReal.push_back(0.0);
            ftImag.push_back(0.0);
            ftAbs.push_back(0.0);
            vec_x.push_back(0.0);
        }
        f_descrete( func , vec_x , dx , number_of_x_sampling );//データの離散化
        for (int k = 0; k < number_of_x_sampling; k++) {
            for (int n = 0; n < number_of_x_sampling; n++) {
                ftReal[k] += vec_x[n] * cos( 2.0 * M_PI / number_of_x_sampling * n * k ) ;
                ftImag[k] += - vec_x[n] * sin( 2.0 * M_PI / number_of_x_sampling * n * k ) ;
            }
            ftAbs[k] = sqrt( ftReal[k] * ftReal[k] + ftImag[k] * ftImag[k] );
        }
        return 0;
    }

    template <typename TInput>
    inline double ftdata(TInput f_x_Data, const double dx, const int number_of_x_sampling){
        for (int i = 0; i < number_of_x_sampling; i++) {//コンストラクタ
            ftReal.push_back(0.0);
            ftImag.push_back(0.0);
            ftAbs.push_back(0.0);
        }
        for (int k = 0; k < number_of_x_sampling; k++) {
            for (int n = 0; n < number_of_x_sampling; n++) {
                ftReal[k] += f_x_Data[n] * cos( 2.0 * M_PI / number_of_x_sampling * n * k ) ;
                ftImag[k] += - f_x_Data[n] * sin( 2.0 * M_PI / number_of_x_sampling * n * k ) ;
            }
            ftAbs[k] = sqrt( ftReal[k] * ftReal[k] + ftImag[k] * ftImag[k] );
        }
        return 0;
    }

};

//---------------------------------------------------------------------------------------------
//    2D 離散フーリエ変換：
//          x-y空間の関数f(x,y)を波数kx(=2π/λ),ky空間の関数F(kx,ky)に変換する。
//
//    ftfメソッド：関数をdxの幅のnumber_of_x_sampling個の離散データに離散化した後に離散フーリエ変換
//    ftfataメソッド：幅dxをもつnumber_of_x_sampling個の離散データを離散フーリエ変換
//    サンプルは、cd ../../program/cpp/test_FT_2D.cpp
//
//---------------------------------------------------------------------------------------------
template <typename TFunc>//関数をdxの幅で離散化する
void f2D_descrete(TFunc func, vector< vector<double> >& vec, const double dx, const double dy,  const int number_of_x_sampling, const int number_of_y_sampling){
    for (int i = 0; i < number_of_x_sampling; i++) {
        for (int j = 0; j < number_of_y_sampling; j++) {
            vec[i][j] = func( i * dx , j * dy );
        }
    }
}

class DFT_2D{//離散フーリエ変換用のクラス
private:
    vector< vector<double> > vec_xy;
public:
    vector< vector<double> > ftReal;
    vector< vector<double> > ftImag;
    vector< vector<double> > ftAbs;

    template <typename TFunc>
    inline double ftf(TFunc func, const double dx, const double dy, const int number_of_x_sampling, const int number_of_y_sampling){
        vector<double> for_constracta1, for_constracta2, for_constracta3, for_constracta4;
        for (int j = 0; j < number_of_y_sampling; j++) {
            for_constracta1.push_back(0.0);
            for_constracta2.push_back(0.0);
            for_constracta3.push_back(0.0);
            for_constracta4.push_back(0.0);
        }
        for (int i = 0; i < number_of_x_sampling; i++) {//コンストラクタ
            ftReal.push_back(for_constracta1);
            ftImag.push_back(for_constracta2);
            ftAbs.push_back(for_constracta3);
            vec_xy.push_back(for_constracta4);
        }

        f2D_descrete( func , vec_xy , dx , dy , number_of_x_sampling , number_of_y_sampling );//データの離散化
        for (int kx = 0; kx < number_of_x_sampling; kx++) {
            for (int ky = 0; ky < number_of_y_sampling; ky++) {
                for (int n = 0; n < number_of_x_sampling; n++) {
                    for (int m = 0; m < number_of_y_sampling; m++) {
                        ftReal[kx][ky] += vec_xy[n][m] * cos( 2.0 * M_PI * n * kx / number_of_x_sampling + 2.0 * M_PI * m * ky / number_of_y_sampling ) ;
                        ftImag[kx][ky] += - vec_xy[n][m] * sin( 2.0 * M_PI * n * kx / number_of_x_sampling + 2.0 * M_PI * m * ky / number_of_y_sampling ) ;
                    }
                }
                ftAbs[kx][ky] = sqrt( ftReal[kx][ky] * ftReal[kx][ky] + ftImag[kx][ky] * ftImag[kx][ky] );
            }
        }
        return 0;
    }

    template <typename TInput>
    inline double ftdata(TInput f_x_Data, const double dx, const double dy, const int number_of_x_sampling, const int number_of_y_sampling){
        vector<double> for_constracta1, for_constracta2, for_constracta3, for_constracta4;
        for (int j = 0; j < number_of_y_sampling; j++) {
            for_constracta1.push_back(0.0);
            for_constracta2.push_back(0.0);
            for_constracta3.push_back(0.0);
        }
        for (int i = 0; i < number_of_x_sampling; i++) {//コンストラクタ
            ftReal.push_back(for_constracta1);
            ftImag.push_back(for_constracta2);
            ftAbs.push_back(for_constracta3);
        }
        for (int kx = 0; kx < number_of_x_sampling; kx++) {
            for (int ky = 0; ky < number_of_y_sampling; ky++) {
                for (int n = 0; n < number_of_x_sampling; n++) {
                    for (int m = 0; m < number_of_y_sampling; m++) {
                        ftReal[kx][ky] += f_x_Data[n][m] * cos( 2.0 * M_PI * n * kx / number_of_x_sampling + 2.0 * M_PI * m * ky / number_of_y_sampling ) ;
                        ftImag[kx][ky] += - f_x_Data[n][m] * sin( 2.0 * M_PI * n * kx / number_of_x_sampling + 2.0 * M_PI * m * ky / number_of_y_sampling ) ;
                    }
                }
                ftAbs[kx][ky] = sqrt( ftReal[kx][ky] * ftReal[kx][ky] + ftImag[kx][ky] * ftImag[kx][ky] );
            }
        }
        return 0;
    }

};

//---------------------------------------------------------------------------------------------
//    高速フーリエ変換：
//          x空間の複素関数f(z)を波数k(=2π/λ)空間の複素関数F(k)に変換する。
//
//    第1引数：fx_data  複素数の座標空間のインプットデータ配列
//    第2引数：Fk_data  複素数の波数空間のアウトプットデータ配列
//
//    サンプルは、cd ../../program/cpp/test_FFT.cpp
//
//---------------------------------------------------------------------------------------------

// inline void FFT(const vector<complex<double>>& fx_data, vector<complex<double>>& Fk_Data){
//     if( ( fx_data.size()&( fx_data.size() - 1 ) ) ) { cout << "Error! The number of input data must be 2^n. " << endl; return;};
//     int n_level;
//     auto& i = n_level;
//     for( i=0; i<64; ++i ){
//         if( fx_data.size()>>i == 1) break;
//     }
//     vector<int> ids;
//     ids.reserve( fx_data.size() );
//     ids.push_back(0);
//     ids.push_back(1);
//     for( int i=0; i < n_level-1; ++i ){
//         auto sz = ids.size();
//         for_each( ids.begin(), ids.end(), [](int& x){ x*=2; } );
//         ids.insert( ids.end(), ids.begin(), ids.end() );
//         auto it = ids.begin();
//         std::advance( it, sz );
//         for_each( it, ids.end(), [](int&x){ x+=1; } );
//     }
//
//     Fk_Data.resize( fx_data.size() );
//     for( int i = 0; i < fx_data.size(); ++i )
//         Fk_Data[i] = fx_data[ids[i]];
//     unsigned int po2 = 1;
//     for( int i_level = 1; i_level <= n_level; ++i_level ){
//         po2 <<= 1;//2倍する
//         const int po2m = po2 >> 1;//1/2倍する
//         auto w = exp( std::complex<double>( 0.0 , 2.0 * M_PI / (double) po2 ) );
//         auto ws = complex<double>( 1.0 , 0.0 );
//         for( int k=0; k<po2m; ++k ){
//             for( int j=0; j < fx_data.size(); j+=po2 ){
//                 auto wfb = ws*Fk_Data[j+k+po2m];
//                 Fk_Data[j+k+po2m] = Fk_Data[j+k] - wfb;
//                 Fk_Data[j+k] = Fk_Data[j+k] + wfb;
//             }ws *= w;
//         }
//     }


}


//------------------------------------------------------------------------------------------------------------
// 　  Brent method
//
//    初期値として極小値を囲い込む3点a,b,cを初期値として選び、その三点を通る放物線で補完する。放物線の最小値を返すxを用いて
//    f(x)=dとして、f(a)f(b)f(c)f(d)の小さいほうから3つを選んで、a,b,cと定義し直して繰り返す。
//
//    brent(ax, bx, cx, func, *xmin, tol): Find a minmum value of a function f(x), and put xmin in adress *xmin.
//    ---> parameter ; func : A function to find a minmum of.
//                     ax,bx,cx : initial guess.
//                     *xmin : a adress where we put the xmin in.
//                     epsilon : (optional) Tolerance for termination.
//------------------------------------------------------------------------------------------------------------

template<typename Tin>
inline Tin sign_for_brent(Tin a, Tin b){
    if ( b >= 0.0 ) return fabs(a);
    else return - fabs(a);
}

template<typename Tin>
inline void shift_for_brent(Tin& a, Tin& b, Tin& c, Tin& d){
    a = b; b = c; c = d;
}

template<typename TFunc , typename Tin>
inline void brent(Tin ax, Tin bx, Tin cx, TFunc func, Tin& xmin, double tol = 1.0e-8){
    constexpr int max_iteration_brent = 200;/* 反復回数の上限 */
    constexpr double cgold = 0.3819660;	/* 黄金分割比 */
    constexpr double zeps_brent = 1.0e-10; /* 極小がちょうど x=0 にあるときは相対精度 tolの代わりにこれを絶対精度とする */

  	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
  	double e = 0.0;						/* 前々回の更新量 */

  	if( ax < cx ){ a = ax ; b = cx; }
  	else { a = cx ; b = ax;}
  	x = w = v = bx;						/* 初期化 */
    fw = fv = fx = func(x);
  	for( int iter = 1 ; iter <= max_iteration_brent; iter++){
    		xm = 0.5 * ( a + b );
    		tol1 = tol * fabs(x) + zeps_brent;
    		tol2 = 2.0 * tol1;
    		if( fabs( x - xm ) <= ( tol2 - 0.5 * ( b - a ) ) ){
      			xmin = x;					/* 最良の値を返す */
    		}
    		if( fabs(e) > tol1 ){
      			r = ( x - w ) * ( fx - fv );
      			q = ( x - v ) * ( fx - fw );
      			p = ( x - v ) * q - ( x - w ) * r;
      			q = 2.0 * (q - r);
      			if(q > 0.0)	p = -p;
      			q = fabs(q);
      			etemp = e;
      			e = d;
      			if(fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
                /* 放物線補間の適否の検査 */
      			{
      				e = (x >= xm)? a - x: b - x;
      				d = cgold * e;			/* 放物線補間は不適。大きい方の区間を */
      			}
      			else{d = p / q;u = x + d;
      				if(u - a < tol2 || b - u < tol2)	d = sign_for_brent(tol1, xm - x);
    			}
    		}

    		else{ e = (x >= xm)? a - x: b - x; d = cgold * e; } //黄金分割

    		u = x + ( ( fabs(d) >= tol1 ) ? d : sign_for_brent( tol1, d ) );
    		fu = func(u);
            					/* 主ループでの関数値評価はここだけ */
    		if(fu <= fx){
      			if(u >= x)	a = x;
      			else		b = x;
      			shift_for_brent( v , w , x , u );
      			shift_for_brent( fv , fw , fx , fu );
    		}
    		else{
      			if( u < x )	a = u;
      			else		b = u;
      			if( fu <= fw || w == x){
        				v = w ; w = u ; fv = fw ; fw = fu;
      			}
      			else if(fu <= fv || v == x || v == w){
        				v = u; fv = fu;
      			}
  		  }
  	}
  	xmin = x;
}


// ---------------------------------------------------------------------------------
//  グラフ描画用のクラス
//  plot( x , y , length )---> x,y :1d-array
//                            length :number of array
//  plotf( func , x_min , x_max , dx )----> func : 1-variable function
//                                        x_min, x_max : range of graph
//                                        dx : (optional) stride of x
//  --example--
//  int main(){ Plot plt; int n = 10; x[n] = {...}; y[n] = {...};
//     plt.plot(x,y,n); plt.xlabel("temp"); plt.show();
//  }
// ---------------------------------------------------------------------------------
// class Plot{
// public:
//     template<typename T>
//     inline void plot(T x[], T y[], int length){
//         ofstream fs;
//         fs.open("graph_nonapy.txt");
//         for (int i = 0; i < length; i++) {
//             fs << x[i] << " " << y[i] << endl;
//         }
//         fs.close();
//     }
//     template<typename TFunc, typename T>
//     inline void plotf(TFunc func, T x_min , T x_max ,T dx = 1.0e-1){
//         ofstream fs;
//         fs.open("graph_nonapy.txt");
//         for (int i = int(x_min/dx); i <= int(x_max/dx); i++) {
//             fs << i * dx << " " << func( i * dx ) << endl;
//         }
//         fs.close();
//     }
//     inline void xlabel(const char* ch1){//x軸のラベルをつけるメソッド
//         strcat(x_label,ch1);
//     }
    // inline void ylabel(const char* ch2){//y軸のラベルをつけるメソッド
    //     strcat(y_label,ch2);
    // }
    // inline void legend(const char* ch3){//凡例をつけるメソッド
    //     strcat(legend_out,ch3);
    // }
    // inline void title(const char* ch4){//タイトルをつけるメソッド
    //     strcpy(title_label,ch4);
    // }
    // inline void title2(const char* ch5){//タイトルをつけるメソッド2
    //     strcpy(title_label2,ch5);
    // }
//     void show(){
//         FILE *pp;
//         pp = popen("python3","w");
//         fprintf(pp, "import pandas as pd \n");
//         fprintf(pp, "import matplotlib.pyplot as plt \n");
//         fprintf(pp, "plt.style.use('ggplot') \n");
//         fprintf(pp, "data=pd.read_csv('graph_nonapy.txt',header=None, delim_whitespace=True) \n");
//         fprintf(pp, "data.plot(x=data.columns[0],y=[data.columns[1]],title='");fprintf(pp,"%s", title_label);fprintf(pp,"')\n");
//         fprintf(pp, "plt.xlabel(' ");fprintf(pp,"%s",x_label);fprintf(pp," ') \n");
//         fprintf(pp, "plt.ylabel(' ");fprintf(pp,"%s",y_label);fprintf(pp," ') \n");
//         fprintf(pp, "plt.legend([' ");fprintf(pp,"%s",legend_out);fprintf(pp," ']) \n");
//         fprintf(pp, "plt.minorticks_on() \n");
//         fprintf(pp, "plt.show() \n");
//         pclose(pp);
//     }
// private:
//     char x_label[100];
//     char y_label[100];
//     char legend_out[100];
//     char title_label[100];
//     char title_label2[100];
// };
//
// };
