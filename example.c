#include "arbcmath.h"

int main()
{
    double complex z;

    z = ac_erfc(20.0 + 10.0*I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));

    z = ac_hyp2f1(1.3 - 2.1*I, 2.5 + 0.6*I, 2.0, 3.0 + 4.0*I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));

    z = ac_hyperu(0.7 - 0.3 * I, 0.7 - 0.2 * I, 10.0 - 100.0 * I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));

    z = ac_hyp1f1(3.14, 2.78, 2015.1130*I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));
}
