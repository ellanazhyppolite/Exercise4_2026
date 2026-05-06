#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include "/Users/ella/Documents/EPFL/BA4/PHYSNUM/Exercise4_2026/Exercise4_2026/common/ConfigFile.h"

using namespace std;

const double PI = 3.1415926535897932384626433832795028841971e0;



// -------------------------------------------------------------------------------------------------------------------------------------

std::valarray<double> k1(2), k2(2), k3(2), k4(2); // preallocate valarrays for the RK4 step parts
std::valarray<double> state(2), trial1(2), halfStep(2), trial2(2), derivativeBuffer(2); // preallocate valarrays for the different steps and dy/dt vector


//int timeScheme = 0;              // 0 = fixed RK4, 1 = adaptive RK4
//int sampling = 1;
//double x = 0.0;
const double xEnd = 1.0;
//double dx = 1.0/64.0;
//double tolerance = 1e-6;
//double d = 100.0;

//double R = 0.0;
//double V0 = 0.0;

int stepsSinceLastOutput = 0;
std::ofstream outputFile;

std::valarray<double> dydx(const std::valarray<double>& state_, double x, double R){
    
    // derivative byffer = dy/dx, suit la formule du latex
    derivativeBuffer[0] = - R * state_[0];
    derivativeBuffer[1] = R + state_[1] / (1 - x);

    return derivativeBuffer;
}

void RungeKutta(const std::valarray<double>& state_, double x, double dx_, double R){
    k1 = dx_*dydx(state_, x, R);
    k2 = dx_*dydx(state_ + 0.5 * k1, x, R);
    k3 = dx_*dydx(state_ + 0.5 * k2, x, R);
    k4 = dx_*dydx(state_ + k3, x, R);
}

double norm(const std::valarray<double>& v){
    double n(0.0);
    for(size_t i = 0; i < v.size() ; ++i){
        n += pow(v[i],2);
    }
    return pow(n,0.5);
}


void step(double x, double dx, double R){

  // en dessous si on fait adptatif
  /* 
    if(timeScheme == 1){
        double d(100.0);
        double tolerance(1e-6);

        do{
            RungeKutta(state, dx, R);
            trial1 = state + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);    //eq 2.141 notes de cours

            RungeKutta(state, dx/2.0, R);
            halfStep = state + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
            RungeKutta(halfStep, dx/2.0, R);
            trial2 = halfStep + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);

            d = norm(trial1 - trial2);

            if(d > tolerance){
                dx = 0.9 * dx * pow(tolerance/d, 1.0/5.0); 
            }
        }while(d > tolerance);

        state = trial1;
        x += dx;

        dx = dx * pow(tolerance/d, 1.0/5.0);

    }
        */

  //  if(timeScheme == 0){
        RungeKutta(state, dx, x, R);
        state += 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
        x += dx;

  //  }
}


double F_alpha(double alpha, double dx, double V0, double R){
    state[0] = V0;
    state[1] = alpha;
    double x(0.0);

    while (x < xEnd) {
        step(x, dx, R);
    } 
    
    return state[1]; // = state[1] - contrainte_initiale (=0)
}

// en dessous printOut
/*
 void printOut(bool force){
        if (force || stepsSinceLastOutput >= sampling) {
            outputFile << x << " ";


        outputFile << state[0] << state[1] << " ";

        outputFile << endl;
        stepsSinceLastOutput = 0;
        
        else {
        stepsSinceLastOutput++;
        }

    }

*/

// -------------------------------------------------------------------------------------------------------------------------------------
   

// Resolution d'un systeme d'equations lineaires par elimination de Gauss-Jordan
// (tridiagonal system: diag, lower, upper, rhs all of consistent sizes)$``

template<class T>
vector<T> solve(const vector<T>& diag,
                const vector<T>& lower,
                const vector<T>& upper,
                const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < (int)diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i]  -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i = (int)diag.size() - 2; i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

// TODO: Implement the relative permittivity epsilon_r(r).
//       Should allow for a trivial test case (trivial=true) 
double epsilon_r(double r, double b, double R, bool trivial)
{
    if(not trivial){
        if (0 <= r and r < b){
            return 1.0;
        }
        else{
            return 3 + 6 * (r - b)/(R - b);
        }
    }
    else{
        return 1.0;
    }
}

// TODO: Implement the normalised free charge density rho_lib(r) / epsilon_0.
//       Should allow for a trivial test case (trivial=true) 
double rho_lib(double r, double b, double a0, bool trivial)
{
    if (not trivial){
        if(0 <= r and r <= b){
            return a0 * sin(PI*r/b);
        }
        else{
            return 0.0;
        }
    }
    else{
        return 1;
    }
}

int main(int argc, char* argv[])
{
    // USAGE: ./engine [configuration-file] [<key=value> ...]

    string inputPath = "trivial.in";
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for (int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    // Physical parameters
    const double b   = configFile.get<double>("b");   // Inner radius [m]
    const double R   = configFile.get<double>("R");   // Outer radius [m]
    const double V0  = configFile.get<double>("V0");  // Boundary potential at r=R [V]
    const double a0  = configFile.get<double>("a0");  // Free charge density scale [V/m^2]
    const bool trivial = configFile.get<bool>("trivial"); // true: uniform test case
    const double dx = configFile.get<double>("dx");

    // Discretisation
    const int N1 = configFile.get<int>("N1"); // Intervals in [0, b]
    const int N2 = configFile.get<int>("N2"); // Intervals in [b, R]

    // Output file prefix
    const string output = configFile.get<string>("output");


    // ---------------------------------------------------------------
    // Build grid
    // ---------------------------------------------------------------
    const int ninters = N1 + N2;         // Total number of intervals
    const int npoints = ninters + 1;     // Total number of grid points
    const double h1 = b / N1;            // Step size in inner region
    const double h2 = (R - b) / N2;      // Step size in outer region

    vector<double> r(npoints);

    // TODO: fill r[0..N1] and r[N1..npoints-1]
    for(size_t i = 0; i < npoints ; ++i){
        if(i <= N1){
            r[i] = i * h1;
        }
        else{
            r[i] = b + (i - N1) * h2;
        }
    }

    vector<double> h(ninters);           // Interval widths
    vector<double> midPoint(ninters);    // Midpoints of each interval

    // TODO: fill h[i]  and  midPoint[i]
    for (size_t i = 0; i < ninters; ++i) {
        h[i] = r[i+1] - r[i];
        midPoint[i] = 0.5 * (r[i] + r[i+1]);
    }


    // ---------------------------------------------------------------
    // Assemble the tridiagonal system  A * phi = rhs
    // ---------------------------------------------------------------
    vector<double> diag(npoints, 0.0);   // Main diagonal
    vector<double> lower(ninters, 0.0);  // Sub-diagonal  (lower[i] links row i+1 to col i)
    vector<double> upper(ninters, 0.0);  // Super-diagonal (upper[i] links row i to col i+1)
    vector<double> rhs(npoints, 0.0);    // Right-hand side

    for (int k = 0; k < ninters; ++k) {
        // TODO: compute alpha_k and beta_k
        double alpha_k(0.0);
        double beta_k(0.0);

        alpha_k = r[k]*epsilon_r(r[k], b, R, trivial) / (h[k]*h[k]);
        beta_k = r[k+1]*epsilon_r(r[k+1], b, R, trivial) / (h[k]*h[k]);

        //       then add their contributions to diag, lower, upper, and rhs
        diag[k] -= alpha_k + beta_k;
        if( k > 0){ lower[k-1] += alpha_k;}
        upper[k] += beta_k;

        rhs[k] -= midPoint[k] * rho_lib(midPoint[k], b, a0, trivial); //le /e_0 déjà dans rho_lib

    }

    // TODO: enforce the Dirichlet BC at r = R
   
   // diag[ninters - 1] = (2 * r[ninters]*epsilon_r(r[ninters], b, R, trivial) + r[ninters - 1]*epsilon_r(r[ninters - 1], b, R, trivial)) / (h2*h2);
   // rhs[ninters - 1] = - midPoint[ninters - 1] * rho_lib(midPoint[ninters - 1], b, a0, trivial) - 2 * midPoint[ninters] * epsilon_r(midPoint[ninters], b, R, trivial) / (h2*h2);

    double beta_last = r[ninters] * epsilon_r(r[ninters], b, R, trivial) / (h2 * h2);
    diag[ninters - 1]  -= beta_last; 
    rhs[ninters - 1]   -= 2.0 * beta_last * V0;
   
    // ---------------------------------------------------------------
    // Solve the linear system
    // ---------------------------------------------------------------
    vector<double> phi = solve(diag, lower, upper, rhs);

    // ---------------------------------------------------------------
    // Compute electric field E_r and displacement D_r (normalised by eps0)
    // ---------------------------------------------------------------
    vector<double> rmid(ninters);
    vector<double> Er(ninters, 0.0);
    vector<double> Dr(ninters, 0.0);
    for (int k = 0; k < ninters; ++k) {
        rmid[k] = midPoint[k];
        
        // TODO: compute Er[k] and Dr[k] 
        Er[k] = - (phi[k+1] - phi[k]) / h[k];
        Dr[k] =  epsilon_r(rmid[k], b, R, trivial) * Er[k];

    }

    // ---------------------------------------------------------------
    // Compute div(D_r)/eps0 and compare to rho_lib/eps0
    // using finite differences on the midpoint values
    // ---------------------------------------------------------------
    vector<double> rmidmid(ninters - 1);
    vector<double> div_Dr(ninters - 1, 0.0);
    vector<double> rho_at_midmid(ninters - 1, 0.0);
    for (int k = 0; k < ninters - 1; ++k) {
        rmidmid[k] = 0.5 * (rmid[k] + rmid[k + 1]);

        // TODO: compute div_Dr[k] and rho_at_midmid[k]
        double rL = rmid[k];
        double rR = rmid[k+1];
        double dr = rR - rL;

        div_Dr[k] = (rR * Dr[k+1] - rL * Dr[k]) / (dr * rmidmid[k]);
        //div_Dr[k] = (Dr[k+1] - Dr[k]) / dr;

        rho_at_midmid[k] = rho_lib(rmidmid[k], b, a0, trivial);
    }

    // ---------------------------------------------------------------
    // Write output files
    // ---------------------------------------------------------------
    {
        // 1. Electric potential: columns  r  phi
        ofstream ofs(output + "_phi.out");
        ofs.precision(15);
        for (int i = 0; i < npoints; ++i)
            ofs << r[i] << " " << phi[i] << "\n";
    }
    {
        // 2. Electric field and displacement: columns  r_mid  E_r  D_r/eps0
        ofstream ofs(output + "_ErDr.out");
        ofs.precision(15);
        for (int k = 0; k < ninters; ++k)
            ofs << rmid[k] << " " << Er[k] << " " << Dr[k] << "\n";
    }
    {
        // 3. Divergence check: columns  r_midmid  div(D_r)/eps0  rho_lib/eps0
        ofstream ofs(output + "_divD_rho.out");
        ofs.precision(15);
        for (int k = 0; k < ninters - 1; ++k)
            ofs << rmidmid[k] << " " << div_Dr[k]
                << " " << rho_at_midmid[k] << "\n";
    }


//  -------------------------------------------------------------------------------------------------------------------------------------

    double alpha0 = 0.0;   // premier guess
    double alpha1 = 1.0;   // deuxième guess
    double tolerance = 1e-6;

    double F0 = F_alpha(alpha0, dx, V0, R);
    double F1 = F_alpha(alpha1, dx, V0, R);

    while (abs(F1) > tolerance) {
        double alpha2 = alpha1 - F1 * (alpha1 - alpha0) / (F1 - F0);
        alpha0 = alpha1; 
        F0 = F1;
        alpha1 = alpha2; 
        F1 = F_alpha(alpha1, dx, V0, R);
    }


    outputFile.open((output + "_tir.out").c_str());
    outputFile.precision(15);

    outputFile << alpha0 << " " << state[0] << " " << state[1] << " " << endl;


    return 0;
}
