# Exam 1
### Shaun Bartschi

## Problem 1
Code implementation is included below in c++.  Please note, that I was unable to derive the correct solution for the second problem of 23, as division by zero seemed to cause several catastrophic errors in my code.
```C++
#include<iostream>
#include<cmath>
#include<vector>

using namespace std;

double f(double x){
        //return sin(x);
        return -sin(x)/x;
}
double df(double x){
        //return cos(x);
        return -(x*cos(x)-sin(x))/(x*x);
}
double ddf(double x){
        //return -sin(x);
        return -((x*x-2)*sin(x)+2*x*cos(x))/(x*x*x);
}

double newtonm(double x0, double tol, int maxiter){
        double error = 10*tol;
        double f0=df(x0);
        double df0 = ddf(x0);
        int cnt = 0;
        double x1;
        int errorCnt = 0;
        double perror = abs(f0/df0);
        while(error > tol && cnt < maxiter){
                if(df0 == 0){
                //      cout<<"invalid denominator"<<endl;
                        x0 = x0 + 0.1;
                        f0 = df(x0);
                        df0 = ddf(x0);
                        continue;
                }
                cnt++;
                x1 = x0-f0/df0;
                error = abs(x1-x0);
                if(error>perror){
                        errorCnt++;
                        if(errorCnt > 5){
                                double round = (df(x1+2*tol)+f(x1-2*tol))/2;
                                if(round <= 10*tol && round >= -10*tol){
                                        return x1;
                                }else{
                                        cout<<"Divergent error at: "<<x1<<endl;
                                        return NULL;
                                }
                        }
                }else{errorCnt=0; perror = error;}
                x0 = x1;
                f0 = df(x0);
                df0 = ddf(x0);
        }
        return x0;
}
int main(){
        double xa = -10;
        double xb = 10;
        int maxiter = 10000;
        double isize = (xb-xa)/1000;
        double tol = pow(10,-8);
        vector<double> xmin; xmin.push_back(xa);
        double fmin = f(xa);
        double oldcrit = xa;
        double testx = xa;
        cout<<"Critical Points:"<<endl;
        for(double i=xa; i<=xb; i+=isize){
                if(df(i)<0.01 && df(i)>-0.01){
                        testx = newtonm(i,tol,maxiter);
                        if(testx > oldcrit+100*tol){
                                double testf = f(testx);
                                if(isinf(testf)){
                                        cout<<"test"<<endl;
                                        testf = (f(testx+tol)+f(testx-tol))/2;
                                }
                                cout<<testx<<"\t"<<testf<<endl;
                                if(ddf(i) >= 0){
                                        double testf = f(newtonm(i,tol,maxiter));
                                        if(testf<fmin){
                                                xmin.clear(); xmin.push_back(testx);
                                                fmin = testf;
                                        }else if (testf == fmin){
                                                xmin.push_back(testx);
                                        }
                                }
                        }
                }
        }
        cout<<"\nMininum Values:"<<endl;
        for(int i=0; i<xmin.size(); i++){
                cout<<xmin[i]<<endl;
        }
        return 0;
}
```

## Problem 2
### 5.
  False, after running a test on my matrix norm functions, a nonsingular matrix failed.
### 7.
  Though I spent some time trying to find a solution to this problem, I was unable to find the solution.
  
## Problem 3  
This is the code in C++:
in main.cpp:
```C++
#include <iostream>
#include "gauss.h"
#include "print.h"
#include <vector>

using namespace std;

int main(){
        vector<vector<double> > A;
        vector<double> x;
        x.push_back(2); x.push_back(1); x.push_back(-1); x.push_back(8);
        A.push_back(x); x.clear();
        x.push_back(-3); x.push_back(-1); x.push_back(2); x.push_back(-11);
        A.push_back(x); x.clear();
        x.push_back(-2); x.push_back(1); x.push_back(2); x.push_back(-3);
        A.push_back(x); x.clear();
        print(A);
        cout<<endl;
        print(gaussJordan(A));
        return 0;

}
```
in gauss.h:
```C++
#ifndef GAUSS_H
#define GAUSS_H
#include <vector>
using namespace std;

double bsub(vector<vector<double> > U, int o){
        int n = U.size();
        vector<double> x; x.resize(n);
        double a=0;
        x[n-1] = U[n-1][n-1]/U[n-2][n-1];
        for(int i = n-2; i>=0; i--){
                for(int j = i+1; j<n; j++){
                        a = a - U[i][j]*x[i];
                        o++;
                }
                x[i] = (U[i][n-1] - a)/U[i][i];
        }
        return o;
}

double gauss(vector<vector<double> > Ab, int o){
        int n = Ab.size();
        double a;
        for(int i = 0; i<n; i++){
                for(int j = i+1; j<n; j++){
                        a = Ab[j][i]/Ab[i][i];
                        for(int k = i; k<n+1; k++){
                                Ab[j][k] -= Ab[i][k]*a;
                                o++;
                        }
                }
        }
        return bsub(Ab, o);
}


double gaussJordan(vector<vector<double> > Ab){
        int m = Ab.size();
        int n = Ab[0].size();
        int o=0;
        double a;
        for(int i = 0; i<m; i++){
                for(int j = 0; j<m; j++){
                        if(i!=j){
                                a = Ab[j][i]/Ab[i][i];
                                for(int k = 0; k<n; k++){
                                        Ab[j][k] -= Ab[i][k]*a;
                                        o++;
                                }
                        }
                }
        }
        return o;
}
#endif
```
Output of the Operation
    gauss: 14
    gauss Jordan: 24
This indicates that Gauss Jordan actually took more that 50% as many computations as gauss method and backsubstitution

## Problem 4
This is the implementation of the method described in problem 17.
```C++
vector<double> telimination(vector<double> a, vector<double> b, vector<double> c, vector<double> x){
        vector<double> o;
        int n = x.size();
        o.resize(n);
        double w = a[0];
        for( int i=1; j<n; i++){
                v[i-1] = c[i-1]/w;
                w = a[i] = b[i]*v[i-1];
                o[i] = (x[i]-b[i]*y[i-1])/w;
        }
        return o;
}
```
