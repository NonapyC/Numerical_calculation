// ---------------------------------------------------------------------------------------------------
// NSLin VERSION 1.0.0
// DIRECTED BY SOSHI NONAKA in 18/06/15
// Copyright (C) 2018 Nonaka Soshi †
// † Master of Science, Department of Physics, Osaka University, Japan.
//
// This is a library for linear algebra  ( library for algebraic computation is nonapy2.hpp ).
//
// E_mail me for inquiry at: sou.235711_at_gmail.com ( exchange: _at_ -> @ )
//---------------------------------------------------------------------------------------------------
//

#include<stdio.h>
#include<math.h>
#include <iostream>
#include <fstream>
using namespace std;

namespace LinNS2{
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// 　　                                   PART2. 線形代数

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------------
//  三重対角行列の　線形連立n次方程式　の解法
//      A_{i,j} * x_{j} = b_{i}という連立方程式をx_{i}解くコード。
//      特に、行列A_{i,j}が三角対角化されている場合はTriDiag関数を用いて連立方程式の解を求めることができる。
//      サンプルコードは ~/../../program/cpp/test_liner.cpp
// ---------------------------------------------------------------------------------

template <typename Matrix, typename Vector>
void TriDiag(Matrix &A, Vector &d, Vector &x, int dim){
    double err = 0.0;
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < i - 1 ; j++) {
            err += A[i][j] + A[j][i];
        }
    }
    if ( err != 0.0 ) {
        cout << "Error : This is not tridiagonal matrix" << endl; //エラー処理
    }else{
        double P[dim], Q[dim];
        for (int i = 0; i < dim; i++) {
            if ( i == 0 ) {
                P[i] = - A[i][i+1] / A[i][i];
                Q[i] = d[i] / A[i][i];
            }else{
                P[i] = - A[i][i+1] / ( A[i][i] + A[i][i-1] * P[i-1] );
                Q[i] = ( d[i] - A[i][i-1] * Q[i-1] ) / ( A[i][i] + A[i][i-1] * P[i-1] );
            }
            x[dim-1] = Q[dim-1];
            for (int i = dim-1; i >= 0; i--) {
                x[i] = P[i] * x[i+1] + Q[i];
            }
        }
    }
}
// ---------------------------------------------------------------------------------
//  任意のn*n行列の LU分解
//      A_{i,j} * x_{j} = b_{i}という連立方程式をx_{i}についてLU分解するコード。
//      LUdecomp関数でLU分解し、そのあとLUSolveで方程式を解く。
//      サンプルコードは ~/../../program/cpp/test_liner.cpp
// ---------------------------------------------------------------------------------
template <typename Matrix>
inline void LUdecomp(Matrix &A, const int dim){
    if(dim <= 0) {cout << "Error : dimention fault" << endl;};

    for(int i = 0; i < dim; ++i){
        for(int j = 0; j <= i; ++j){// L_ijの計算(i >= j)
            double lu = A[i][j];
            for(int k = 0; k < j; ++k){
                lu -= A[i][k] * A[k][j];    // L_ik * U_kj
            }
            A[i][j] = lu;
        }
        for(int j = i+1; j < dim; ++j){// U_ijの計算(i < j)
            double lu = A[i][j];
            for(int k = 0; k < i; ++k){
                lu -= A[i][k] * A[k][j];    // L_ik * U_kj
            }
            A[i][j] = lu / A[i][i];
        }
    }
}
// ---------------------------------------------------------------------------------
//  任意のn*n行列のLU分解を用いた　線形連立n次方程式　の解法
//      サンプルコードは ~/../../program/cpp/test_liner.cpp
// ---------------------------------------------------------------------------------
template <typename Matrix, typename Vector>
inline void LUSolve(const Matrix &A, const Vector &b, Vector &x, const int dim){
    if(dim <= 0) {cout << "Error : dimention fault" << endl;};

    for(int i = 0; i < dim; ++i){//  LY=bからYを計算
        double bly = b[i];
        for(int j = 0; j < i; ++j){
            bly -= A[i][j]*x[j];
        }
        x[i] = bly/A[i][i];
    }
    for(int i = dim-1; i >= 0; --i){//  UX=YからXを計算
        double yux = x[i];
        for(int j = i+1; j < dim; ++j){
            yux -= A[i][j]*x[j];
        }
        x[i] = yux;
    }
}

// ---------------------------------------------------------------------------------
//  ガウスの消去法 (ピボット交換なし)
//      A : n×nの係数項とn×1の定数項(b)を併せたn×(n+1)の行列．n+1列目に解が入る
//      dim ; 行列 Aの行の次元 n.
//      サンプルコードは ~/../../program/cpp/test_liner.cpp
// ---------------------------------------------------------------------------------
template<typename Matrix>
void GaussElimination(Matrix &A, int dim){
    // 前進消去(forward elimination)で A+b の拡大行列を上三角行列にする。
    for(int k = 0; k < dim-1; ++k){
        double akk = A[k][k];
        for(int i = k+1; i < dim; ++i){
            double aik = A[i][k];
            for(int j = k; j < dim+1; ++j){
                A[i][j] = A[i][j] - aik * ( A[k][j] / akk ) ; //i行目のj列成分を0にするよう更新していく。
            }
        }
    }
    // 後退代入(back substitution)
    A[dim-1][dim] = A[dim-1][dim]/A[dim-1][dim-1]; //漸化式の末項
    for(int i = dim-2; i >= 0; --i){
        double ax = 0.0;
        for(int j = i+1; j < dim; ++j){
            ax += A[i][j] * A[j][dim];
        }
        A[i][dim] = ( A[i][dim] - ax ) / A[i][i]; //漸化式を解く。
    }
}

};
