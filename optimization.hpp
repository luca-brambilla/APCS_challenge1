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

namespace optimization{

    typedef std::vector<double> Point;
    typedef std::vector<double> Grad;
    
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
                            const Point &x, const double &alpha0,
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
        Point xk(x);
        Point xtmp;
        Point grad( df(xk) );
        
        //! done twice - find a way to work with while loop
        for (size_t i = 0; i < xk.size(); ++i)
        {
            xtmp[i] = xk[i] - a * grad[i];
        }
        double residual = f(xk) - f(xtmp);
        double vec_norm = norm(grad);

        while ( residual >= sigma * a * vec_norm * vec_norm) 
        {
            std::cout << a << std::endl;
            a /= 2.0;
            for (size_t i = 0; i < xk.size(); ++i)
            {
                xtmp[i] = xk[i] - a * grad[i];
            }
            residual = f(xk) - f(xtmp);
            vec_norm = norm(grad);
        }

        return a;
    }

    // template <class Function, class Dfunction>
    // Point gd(Function const &f, Dfunction const &df,
    inline Point gd(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr)
        /**
        * Minimize using gradient descent with adaptive learning rate
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
        Point xtmp(x0);
        Point grad(x0);

        //! put infinity
        double vec_norm = tols+1.0;
        double residual = tolr+1.0;
        double alpha = alpha0;

        // begin iteration loop
        for (unsigned int k=0; (vec_norm>tols) or (residual>tolr) or (k>maxiter); ++k )
        {
            std::cout << "iteration: " << k << std::endl;
            // update gradient and learning rate
            grad = df(xk);
            alpha = armijo(f,df,xk,alpha);

            // gradient descent element-wise
            for (size_t i=0; i<xk.size(); ++i)
            {
                xtmp[i] = xk[i];
                xk[i] -= alpha * grad[i];
            }

            // update gradient and learning rate
            grad = df(xk);
            alpha = armijo(f,df,xk,alpha0);

            residual = std::abs( f(xk) - f(xtmp) );
            vec_norm = norm(grad);

            optimization::print(xk);
        }

        return xk;
    }

    // template <class Function, class Dfunction>
    // Point optimize(Function const &f, Dfunction const &df,
    inline Point optimize(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr)
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
        Point res = gd(f,df,x0,alpha0,maxiter,tols,tolr);

        return res;
    }


} // namespace optimization

#endif