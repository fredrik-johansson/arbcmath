#include "arbcmath.h"

int main()
{
    double complex z;

    z = ac_erfc(20.0 + 10.0*I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));

    z = ac_hyp2f1(1.3 - 2.1*I, 2.5 + 0.6*I, 2.0, 3.0 + 4.0*I);

    printf("%.15g + %.15g*I\n", creal(z), cimag(z));
}
