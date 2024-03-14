// Luca Brambilla
// 10510718 - 919812

#include <cstddef>
#include <functional>
#include <vector>
#include <string>
#include <cmath>            // exp, log
#include <iostream>         // input output
//#include <GetPot>           // read data from file
//#include "muparser_fun.hpp" // read functions from string

#include "optimization.hpp"

double myfun(std::vector<double> x)
{
    return x[0]*x[1] + 4*x[0]*x[0]*x[0]*x[0] + x[1]*x[1] + 3*x[0];
}
optimization::Point myfun_grad(std::vector<double> x)
{
    std::vector<double> grad(2);
    grad[0] = x[1] + 16.0*x[0]*x[0]*x[0] + 3.0;
    grad[1] = x[0] + 2.0*x[1];
    return grad;
}

optimization::Point gradient(std::vector<optimization::Function> grad, optimization::Point x)
{
    optimization::Point res(x);
    if (grad.size() != x.size())
    //! return not a number
        return res;

    for(unsigned int i=0; i<x.size(); ++i)
    {
        res[i] = grad[i](x);
    }

    return res;
}

int main(int argc, char **argv)
{
    double epsr = 1e-6;
    double epss = 1e-6;
    double alpha0 = 1.0;
    unsigned int maxiter = 100;
    optimization::Point x0(2);
    optimization::Method method = optimization::moment;

    //! remove optimization and do it in main directly
    optimization::Point x_opt( optimization::optimize(myfun, myfun_grad,x0,alpha0,maxiter,epss,epsr,method) );
    optimization::print(x_opt);
}