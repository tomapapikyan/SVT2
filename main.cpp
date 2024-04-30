#include <iostream>
#include <cmath>
#include "inmost.h"


using namespace INMOST;
using namespace std;

void
get_data(int N, Sparse::Matrix& A, Sparse::Vector& b, Sparse::Vector& exact_sol){

    for(int i=0; i<(N-1)*(N-1); i++) {
        A[i][i] = 4;
        exact_sol[i] = cos(5.0*(i%(N-1) + 1)/N)*sin(5.0*(i/(N-1) + 1)/N);
        b[i] = 50*cos(5.0*(i%(N-1) + 1)/N)*sin(5.0*(i/(N-1) + 1)/N)/N/N;

        if(i%(N-1)==0)
            b[i] += sin(5.0*(i/(N-1) + 1)/N);
        else
            A[i][i-1] = -1;
                    if(i%(N-1) == N-2)
            b[i] += cos(5)*sin(5.0*(i/(N-1) + 1)/N);
        else
            A[i][i+1] = -1;

        if(i/(N-1)==N-2)
            b[i] += cos(5.0*(i%(N-1) + 1)/N)*sin(5);
        else
            A[i][i+N-1] = -1;

        if (i/(N-1) > 0)
            A[i][i-N+1] = -1;
    }
}

void
test_num_sol(int N, Sparse::Vector num_sol, Sparse::Vector exact_sol)
{
    double val, l2_norm = 0, c_norm = 0;

    for(int i=0; i<(N-1)*(N-1); i++) {
        val = fabs(num_sol[i] - exact_sol[i]);
        l2_norm += val * val;

        if (val > c_norm) {
            c_norm = val;
        }
    }

    cout << "l2_norm: " << sqrt(l2_norm/((N-1)*(N-1))) << endl;
    cout << "c_norm: " << c_norm<< endl;
}




int main()
{
    int N;
    cin >> N;

    Sparse::Matrix A;
    Sparse::Vector b;
    Sparse::Vector num_sol;
    Sparse::Vector exact_sol;

    A.SetInterval(0, (N-1)*(N-1));
    b.SetInterval(0, (N-1)*(N-1));
    num_sol.SetInterval(0, (N-1)*(N-1));
    exact_sol.SetInterval(0, (N-1)*(N-1));

    get_data(N, A, b, exact_sol);

    Solver S(Solver::INNER_MPTILU2);

    S.SetParameter("absolute_tolerance", "1e-10");
    S.SetParameter("relative_tolerance", "1e-6");

    S.SetMatrix(A);
   bool solved = S.Solve(b, num_sol);

    cout << "num.iters: " << S.Iterations() << endl;
    cout << "prec.time: " << S.PreconditionerTime() << endl;
    cout << "iter.time: " << S.IterationsTime() << endl;

    if(!solved){
        cout << "Linear solver failure!" << endl;
        cout << "Reason: " << S.ReturnReason() << endl;
    }

    test_num_sol(N-1, num_sol, exact_sol);
    return 0;
}

