// Luca Brambilla
// 10510718 - 919812

#include <vector>
#include <string>
#include <cmath>            // exp, log
#include <iostream>         // input output
#include <GetPot>           // read data from file
#include "muparser_fun.hpp" // read functions from string
#include "solver_theta.hpp" // ODE solver with theta method

// GNUPLOT variable must be defined while compiling if 
// gnuplot and gnuplot-iostream are available
#ifdef GNUPLOT
#include "gnuplot-iostream.hpp" // interface with gnuplot
#endif


// function definition alternative to muparser

// double myfun(double t, double y)
// {
//     return -t*std::exp(-y);
// }
// double myfun_d(double t, double y)
// {
//     return t*std::exp(-y);
// }
// double myfun_ex(double t)
// {
//     return std::log(-t*t/2+1);
// }

double my_err(std::vector<double> const &v1, std::vector<double> const &v2)
{
    /*
    computes the error between two vectors in infinity norm
    @param v1       vector 1
    @param v2       vector 2

    @return norm    norm of the error
    */
    if (v1.size() != v2.size())
    {
        std::cout << "vector are of different size" << std::endl;
        return 0.;
    }

    double norm = 0.;
    double tmp = 0.;
    for (size_t i=0; i<v1.size(); ++i)
    {
        //norm += std::abs( v1[i]-v2[i] );
        tmp = std::abs(v1[i]-v2[i]);
        tmp > norm ? norm = tmp : 0;
    }

    return norm;
}


int main(int argc, char **argv)
{

// ! PARSE COMMAND LINE WITH GETPOT
    GetPot indata(argc, argv);
    if (indata.search(2, "--help", "-h"))
    {
        std::cout << "Usage: " << argv[0] << "<-i DataFile>" << std::endl;
        std::cout << "DataFile is a GetPot file. If missing is set to data.pot"
            << std::endl;
        return 0;
    };
    std::string dataFile;
    dataFile = indata.follow("data.pot", "-i");
    std::cout << " Reading data from " << dataFile << std::endl;

// ! READ PARAMETERS FROM GETPOT
    GetPot         ifl(dataFile.c_str());

    const double y0               = ifl("y0", 0.0);
    const double T                = ifl("T", 1.0);
    const double theta            = ifl("theta", 0.5);
    const unsigned int N_INIT     = ifl("N_INIT", 2);
    const unsigned int FACTOR     = ifl("FACTOR", 5);
    const unsigned int NUM_ERR    = ifl("NUM_ERR", 6);

    const std::string fun_str     = ifl("fun_str", "-t*exp(-y)");
    const std::string dfun_str    = ifl("dfun_str", "t*exp(-y)");
    const std::string exfun_str   = ifl("exfun_str", "log(-t*t/2+1)");
    
//! USE MU PARSER FUNCTIONS
    MuparserFun2D myfun(fun_str), myfun_d(dfun_str);
    MuparserFun1D myfun_ex(exfun_str);

//! SOLVE THE CAUCHY PROBLEM
    unsigned int N = N_INIT;
    double h = 1.;
    double t = 0.0;

    if ( (theta < 0) || (theta > 1) )
    {
        std::cout << "theta must be in [0,1]" << std::endl;
        return 1;
    }

    // initialize vector to store data for convergence analysis
    std::vector<double> err(NUM_ERR);
    std::vector<double> vh(NUM_ERR);
    
    // at each iteration a different N is used
    for (unsigned int k = 0; k<NUM_ERR; ++k)
    {
        N *= FACTOR;
        std::vector<double> y_h(N + 1);
        std::vector<double> coor(N + 1);
        std::vector<double> y_ex(N + 1);

        y_h = solver(myfun, myfun_d, y0, T, N, theta);

        h = T / static_cast<double>(N);
        vh[k] = h;


// ! PRINT ODE RESULTS TO FILE

        // writing results with format
        // t_i, y_h(t_i), y(t_i) and launch gnuplot

        std::cout << "Result file: result_N"<< N <<".dat" << std::endl;
        std::string name = "./output/result_N";
        name.append(std::to_string(N)).append(".dat");
        std::ofstream f(name);
        f << "#t_i\t\tu_h\t\t\tu_ex" << std::endl;

        t = 0.0;
        for(unsigned int i = 0; i <= N; ++i)
        {
            f.setf(std::ios::left);
            f.setf(std::ios::showpos);
            f.setf(std::ios::scientific);
            f.setf(std::ios::showpoint);
            f.precision(3);
            f.width(10);

            coor[i] = t;
            y_ex[i] = myfun_ex(t);

            f << coor[i] << "\t" << y_h[i] << "\t" << y_ex[i] << "\n";
            t += h;
        }

        err[k] = my_err(y_h, y_ex);
        
        f.close();

        // ! PLOT SOLUTION OF ODE CHANGING N FROM VARIABLES
#ifdef GNUPLOT
#ifdef GNUTEMP
        // gnuplot iostream
        Gnuplot gp;
        // naming axes and title
        gp << "set xlabel 't'; set ylabel 'y'; " <<
            "set title 'N = " << std::to_string(N) << ", h = "<< std::to_string(vh[k]) << "'" << std::endl;
        // use variables and temp file
        gp << "plot" << gp.file1d(std::tie(coor, y_h)) << "w lp lw 5 title 'u_h',"
            << gp.file1d(std::tie(coor, y_ex)) << "w l lw 5 title 'u_{ex}'"
            << std::endl;
#endif
#endif

    }

// ! PRINT CONVERGENCE OF EACH CASE

    std::cout << "Result file: convergence.dat" << std::endl;
    std::ofstream g("./output/convergence.dat");
    g << "#h\t\t\terr\t\t\ttheo" << std::endl;
    N=N_INIT;
    double theo = 0.0;
    for(unsigned int k = 0; k < NUM_ERR; ++k)
    {
        g.setf(std::ios::left);
        g.setf(std::ios::showpos);
        g.setf(std::ios::scientific);
        g.setf(std::ios::showpoint);
        g.precision(3);
        g.width(10);
        N *= FACTOR;

        if (theta == 0.5){
            // crank-nicolson
            theo = vh[k]*vh[k]/vh[0]/vh[0]*err[0];
        } else {
            // euler or other
            theo = vh[k]/vh[0]*err[0];
        }

        g << vh[k] << "\t" << err[k] << "\t" << theo <<"\n";

        // ! PLOT SOLUTION OF ODE CHANGING N FROM FILE
#ifdef GNUPLOT    
#ifdef GNUPERM
        Gnuplot gp; // gnuplot iostream

        // preprocessing
        //gp << "load plot.plt" << std::endl;

        // naming axes and title
        gp << "set xlabel 't'; set ylabel 'y'; " <<
            "set title 'N = " << std::to_string(N) << ", h = "<< std::to_string(vh[k]) << "'" << std::endl;

        // read data from existing file
        gp << "plot './output/result_N"<< std::to_string(N) <<".dat' u 1:2 w lp lw 5 title 'u_h', "
        << "'./output/result_N"<< std::to_string(N) <<".dat' u 1:3 w l lw 5 title 'u_{ex}'"
        << std::endl;

        #ifdef PLOTSAVE
        // postprocessing and saving plot
        gp << "set size ratio 0.66; set terminal pdf color font ',25' size 9,6;"
            << "set output './figures/plot_N" << std::to_string(N) << ".pdf';"
            << "replot; set term x11 " << std::endl;
        #endif

#endif
#endif

    }
    g.close();


    // ! PLOT CONVERGENCE ANALYSIS
#ifdef GNUPLOT    
    Gnuplot gp;

    // preprocessing
    gp << "set logscale xy" << std::endl;
    gp << "set key top left" << std::endl;
    gp << "set xlabel 'h'; set ylabel 'err'; "
        << "set title 'theta = " << std::to_string(theta) <<"'" << std::endl;
    gp << "plot './output/convergence.dat' u 1:2 w lp lw 5 title 'err', ";

    // change legend label according to theta
    if (theta == 0.5) {
        gp << "'./output/convergence.dat' u 1:3 w l lw 5 title 'h^2'" << std::endl;
    } else {
        gp << "'./output/convergence.dat' u 1:3 w l lw 5 title 'h'" << std::endl;
    }


    #ifdef PLOTSAVE
    // postprocessing and saving plot
    gp << "set size ratio 0.66; set terminal pdf color font ',25' size 9,6;"
    << "set output './figures/conv_theta" << std::to_string(theta) << ".pdf';"
    << "replot; set term x11 " << std::endl;
    #endif

#endif

    return 0;
}


