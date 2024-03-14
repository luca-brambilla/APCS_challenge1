// Luca Brambilla
// 10510718 - 919812

#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

#include <iostream>
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

//! typedef vector, can also be a different container?
//! initialize vectors with dimensions instead of copy constructor?
//! define direction-wise derivatives to optimize computation

namespace optimization{

    typedef std::vector<double> Point;
    typedef std::vector<double> Grad;
    typedef std::function<double(Point)> Function;
    typedef std::function<Point(Point)> Gradient;

    enum Method {gd, moment};

    struct parameters{
        Point x0;
        double alpha0 = 1.0;
        unsigned int maxiter = 100;
        double tols = 1e-6;
        double tolr = 1e-6;
        double sigma = 0.2;
        double eta = 0.9;
    };

    inline void print(Point const &x)
    {
        for (auto it = x.cbegin(); it != x.cend(); ++it)
        {
            std::cout << *it << " " << std::endl;
        }
        std::cout << std::endl;
    }

    //! why inline???
    inline double norm(Point const &vec)
        /**
        * Compute the euclidean norm of a vector
        * 
        * @param vec      vector
        * 
        * @return res     euclidean norm of the vector
        **/
    {
        double res = 0.0;
        for (auto it = vec.cbegin(); it != vec.cend(); ++it)
        {
            res += *it * *it;
        }

        res = std::sqrt(res);

        return res;
    }

    //template <class Function, class Dfunction>
    //double armijo(Function const &f, Dfunction const &df,
    inline double armijo(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &xk, const double &alpha0,
                            const double &sigma = 0.25)
        /**
        * Armijo rule for adaptive learning rate
        * 
        * @param f        function to be optimized
        * @param df       derivative of the function to be optimized
        * @param x0       initial value
        * @param alpha0   initial learning rate
        * @param sigma    coefficient for stopping criterion
        * 
        * @return a       double containing the optimal learning rate
        **/
    {
        double a = alpha0;
        Point xtmp(xk);
        Point grad( df(xk) );
        
        //! done twice - find a way to work with while loop
        for (size_t i = 0; i < xk.size(); ++i)
        {
            xtmp[i] = xk[i] - a * grad[i];
        }
        double residual = f(xk) - f(xtmp);
        double vec_norm = norm(grad);

        // check for stopping rule
        while ( residual <= sigma * a * vec_norm * vec_norm) 
        {
            // update learning rate
            a /= 2.0;
            
            // check for new residual
            for (size_t i = 0; i < xk.size(); ++i)
            {
                xtmp[i] = xk[i] - a * grad[i];
            }
            residual = f(xk) - f(xtmp);
        }

        return a;
    }

    // template <class Function, class Dfunction>
    // Point gd(Function const &f, Dfunction const &df,
    inline Point gradient_descent(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr)
        /**
        * Minimize using the gradient descent method with adaptive learning rate
        * 
        * @param f        function to be optimized
        * @param df       derivative of the function to be optimized
        * @param x0       initial value
        * @param alpha0    learning rate
        * @param maxiter  maximum number of iterations
        * @param tols     tolerance for the step length
        * @param tolr     tolerance for the residual
        * 
        * @return res     double containing the optimization point
        **/
    {
        // initialize points with correct dimensions
        Point xk(x0);
        Point xtmp(x0.size());
        Point grad(x0.size());

        //! put infinity
        double vec_norm = tols+1.0;
        double residual = tolr+1.0;
        double alpha = alpha0;

        // begin iteration loop
        for (unsigned int k=0; (vec_norm>tols) and (residual>tolr) and (k<maxiter); ++k )
        {
            // std::cout << "iteration: " << k << std::endl;
            // update gradient and learning rate
            grad = df(xk);
            alpha = armijo(f,df,xk,alpha0);

            // gradient descent element-wise
            for (size_t i=0; i<xk.size(); ++i)
            {
                xtmp[i] = xk[i];
                xk[i] -= alpha * grad[i];
            }

            residual = std::abs( f(xk) - f(xtmp) );
            vec_norm = norm(grad);
        }

        return xk;
    }

// template <class Function, class Dfunction>
    // Point gd(Function const &f, Dfunction const &df,
    inline Point momentum(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr,
                            const double eta=0.9)
        /**
        * Minimize using the momentum method with adaptive learning rate
        * 
        * @param f        function to be optimized
        * @param df       derivative of the function to be optimized
        * @param x0       initial value
        * @param alpha0   learning rate
        * @param maxiter  maximum number of iterations
        * @param tols     tolerance for the step length
        * @param tolr     tolerance for the residual
        * @param eta      memory parameter
        * 
        * @return res     double containing the optimization point
        **/
    {
        // initialize points with correct dimensions
        Point xk(x0);
        Point xtmp(x0);
        Point grad(x0.size());
        Point d(x0.size());

        //! put infinity
        double vec_norm = tols+1.0;
        double residual = tolr+1.0;
        double alpha = alpha0;

        // begin iteration loop
        for (unsigned int k=0; (vec_norm>tols) and (residual>tolr) and (k<maxiter); ++k )
        {
            std::cout << "iteration: " << k << ", maxiter: " << maxiter << std::endl;
            // update gradient
            grad = df(xk);

            // momentum method element-wise
            for (size_t i=0; i<xk.size(); ++i)
            {
                // update step
                d[i] = eta*d[i] - alpha*grad[i];
                // save old point
                xtmp[i] = xk[i];
                // find new point
                xk[i] += d[i];
            }
            print(xk);

            // learning rate
            alpha = armijo(f,df,xk,alpha0);

            residual = std::abs( f(xk) - f(xtmp) );
            vec_norm = norm(grad);
        }

        return xk;
    }

    // template <class Function, class Dfunction>
    // Point optimize(Function const &f, Dfunction const &df,
    inline Point optimize(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr,
                            Method m)
        /**
        * Minimize using a chosen method among ... with adaptive learning rate
        * 
        * @param f        function to be optimized
        * @param df       derivative of the function to be optimized
        * @param x0       initial value
        * @param alpha0   inital learning rate
        * @param maxiter  maximum number of iterations
        * @param tols     tolerance for the step length
        * @param tolr     tolerance for the residual
        * 
        * @return res     double containing the optimization point
        **/
    {
        Point res(x0.size());
        switch(m)
        {
            case gd : 
                res = gradient_descent(f,df,x0,alpha0,maxiter,tols,tolr);
                std::cout << "Gradient descent method" << std::endl;
                break;
            case moment : 
                res = momentum(f,df,x0,alpha0,maxiter,tols,tolr);
                std::cout << "Momentum method" << std::endl;
                break;
        }
        return res;
    }


} // namespace optimization

#endif