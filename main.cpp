#include <iostream>
#include <math.h>
#include "rkf45.h"
#include "quanc8.h"
#include "zeroin.h"
#include "rkf45.cpp"
#include "zeroin.cpp"
#include "quanc8.cpp"

double E;
double y[2];
double dy[2];

// подынтигральная функция
double integral(double t){
    double i = (t*t+1)*(3*t*t+4);
    return 1/sqrt(i);
}

double equation(double x){
    return x*x - cos(x);
}

void system(double t, double *y, double *dy){
    dy[0] = y[1];
    dy[1] = -(1 + E * y[0] * y[0]) * y[0];
}

double abserr = 0, relerr = 0;
double result, errest, flag;
int nofun;
double tol = 0.00001;
double x, A, B = 0;

int n = 2;
double h = 0.4;
double t, tout = 0;
double re = 1.0e-5, abser = 0;
int iflag = 1;
double work[15];
int iwork[30];

void init45() {
    iflag = 1;
    t = 0;
}

int main() {
    quanc8(&integral, 0, 1, abser, relerr, &result, &errest, &nofun, &flag);
    E = result * 0.497286;
    printf("E = %.8f\n", E);
    x = Zeroin(equation, 0, 1, tol);
    A = 1.213399 * x;
    printf("A = %.8f\n", A);
    printf("\n");
    printf(" t |  u(t)  \n");
    printf("----------------\n");
    y[0] = A;
    y[1] = B;

    for (int i = 1; i <= 16.1 / h; i++) {
        init45();
        tout = h * i;
        t = tout - h;
        RKF45(system, n, y, &t, &tout, &re, &abser, &iflag, work, iwork);
        printf("%3.1f| ", tout);
        printf("%.8f  ", y[0]);
        printf("\n");
    }

    return 0;
}

