// Luca Brambilla
// 10510718 - 919812


#include <functional>
#include <vector>

double nonlin_fun(std::function<double(double,double)> const &f, double x, double t_n_1, double h, double u_n, double theta)
{
    /*
    Non linear function - part of the Newton scheme to solve ODE

    @param f        function of the Cauchy problem
    @param x        solution approximation at next time step u_{n+1}
    @param t_n_1    next time step t_{n+1}
    @param h        time discretization
    @param u_n      solution approximation from previous step u_{n}
    @param theta    theta for theta-method

    @return 
    */
    
    //! t_n_1-h is instant t_{n}
    return x - u_n - ( (1-theta)*f(t_n_1, x) + theta*f(t_n_1-h, u_n) )*h;
}

double nonlin_fun_d(std::function<double(double,double)> const &df, double x, double t_n_1, double h, double theta)
{
    /*
    Derivative of non linear function - part of the Newton scheme to solve ODE

    @param df       derivative function of the Cauchy problem
    @param x        solution approximation at next time step u_{n+1}
    @param t_n_1    next time step t_{n+1}
    @param h        time discretization
    @param theta    theta for theta-method

    @return 
    */

    return 1 - (1-theta)*df(t_n_1, x)*h;
}

template <class Function, class Dfunction, class NLfunction, class DNLfunction>
std::tuple<double, bool> Newton(Function const &f, Dfunction const & df,
                                NLfunction const &nlf, DNLfunction const &dnlf,
                                double a, double t1, double h, double u_n, double theta,
                                double tol = 1.e-4, double tola = 1.e-10, unsigned int maxIt = 150)
{
    /*
    Newton scheme to find iteratively a zero of a non linear equation 

    @param f        function of the Cauchy problem to be passed to nlf
    @param df       derivative of the function of the Cauchy problem to be passed to dnlf
    @param nlf      nonlinear function to be solved with Newton method
    @param dnlf     derivative of the nonlinear function to be solved with Newton method
    @param a        initial guess for u_{n+1}
    @param t1       time instant t_{n+1} for the next solution u_{n+1}
    @param h        discretization step
    @param u_n      previous solution
    @param theta    value to be passed for the theta-method
    @param tol      tolerance
    @paral tola     absolute tolerance
    @param maxIt    maximum number of iteration

    @return         tuple containing the u_{n+1} and a boolean indicating if maxIt was reached
    */


    double       ya = nlf(f, a, t1, h, u_n, theta);
    double       resid = std::abs(ya);
    unsigned int iter{0u};
    double       check = tol * resid + tola;
    bool         goOn = resid > check;
    while(goOn && iter < maxIt)
    {
        ++iter;
        a += - ya/dnlf(df, a, t1, h, theta);
        ya = nlf(f, a, t1, h, u_n, theta);
        resid = std::abs(ya);
        goOn = resid > check;
    }

    // std::cout << iter << std::endl;

    if (iter > maxIt)
    {
        std::cout << iter << std::endl;
    }

    return std::make_tuple(a, (iter < maxIt));
}


template <class Function, class Dfunction>
std::vector<double> solver(Function const &f, Dfunction const &df,
                           const double &y0, const double &T, const unsigned int &N, const double &theta)
{
    /*
    ODE solver using the theta-method - exploits Newton scheme for implicit cases

    @param f        function of the Cauchy problem
    @param df       derivative of the function of the Cauchy problem
    @param y0       initial value
    @param T        right extremum of the integration integral
    @param N        number of time steps
    @param theta    value for the theta-method

    @return res     vector containing the time history of the solution
    */


    // time step
    double h = T / static_cast<double>(N);

    // store results in vector
    std::vector<double> res;
    res.reserve(N+1);
    // store y0
    res.emplace_back(y0);

    // zeros of nonlinear function at each time step
    std::tuple<double,bool> u_n_1;

    // instant t_{n+1}
    double t_n_1 = 0.0;

    for (unsigned int i=0; i<N; ++i)
    {
        // step forward in time
        t_n_1 += h;
        // solution of nonlinear equation
        // set as initial guess the approximation from previous step u_{n}
        u_n_1 = Newton(f, df, nonlin_fun, nonlin_fun_d, res[i], t_n_1, h, res[i], theta);
        // store u_{n+1}
        res.emplace_back(std::get<0>(u_n_1));

    }

    return res;
}