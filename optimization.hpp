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
//! check pass by reference and how to pass vectors

namespace optimization{

    typedef std::vector<double> Point;
    typedef std::vector<double> Grad;
    typedef std::function<double(Point)> Function;
    typedef std::function<Point(Point)> Gradient;

    enum Method {gd, moment, ad};
    enum Decay {arm, exp, inv};

    struct Parameters{
        Point x0;
        double alpha0 = 1.0;
        unsigned int maxiter = 100;
        double tols = 1e-6;
        double tolr = 1e-6;
        double coeff = 0.2; // for adaptive learning rate
        double eta = 0.9;   // for optimizazion
                            // beta1 for adam
        double beta2 = 0.999; // for adam
        double eps = 1.0e-8;  // for adam
        std::function<double(Point)> f;
        std::function<Point(Point)> df;
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

    inline double armijo(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &xk, const double &alpha0,
                            const double &sigma = 0.25);
    inline double inverse(size_t const &k,
                            const double &alpha0,
                            const double &mu = 0.2);
    inline double exponential(size_t const &k,
                            const double &alpha0,
                            const double &mu = 0.2);                            

    //template<enum Decay>
    double compute_decay(Decay const &dec, const size_t &k,
                            const Point &xk, const double &alpha0,
                            const double &coeff = 0.25,
                            std::function<double(Point)> f = [](Point x){return x[0];},
                            std::function<Point(Point)> df = [](Point x){return x;})
    {
        double alpha = 1.;
        if (dec == Decay::arm)
        {
            alpha = armijo(f, df, xk, alpha0, coeff);
        }
        else if (dec == Decay::inv)
        {
            alpha = inverse(k, alpha0, coeff);
        }
        else if (dec == Decay::exp) {
            alpha = exponential(k, alpha0, coeff);
        }

        std::cout<< alpha << std::endl;
        return alpha;
    }

    //template <class Function, class Dfunction>
    //double armijo(Function const &f, Dfunction const &df,
    inline double armijo(std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &xk, const double &alpha0,
                            const double &sigma)
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

    inline double exponential(size_t const &k,
                            const double &alpha0,
                            const double &mu)
        /**
        * Exponential decay rule for adaptive learning rate
        * 
        * @param k        iteration
        * @param alpha0   initial learning rate
        * @param mu       exponent
        * 
        * @return a       double containing the optimal learning rate
        **/
    {
        return alpha0*std::exp(-mu*k);
    }

    inline double inverse(size_t const &k,
                            const double &alpha0,
                            const double &mu)
        /**
        * Exponential decay rule for adaptive learning rate
        * 
        * @param k        iteration
        * @param alpha0   initial learning rate
        * @param mu       denominator coefficient
        * 
        * @return a       double containing the optimal learning rate
        **/
    {
        return alpha0/(1-mu*k);
    }

    // template <class Function, class Dfunction>
    // Point gd(Function const &f, Dfunction const &df,
    inline Point gradient_descent( const Decay &dec, 
                            std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr, const double& coeff)
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
            alpha = compute_decay(dec,k,xk,alpha0,coeff,f,df);

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
    inline Point momentum( const Decay &dec,
                            std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr, const double &coeff=0.9,
                            const double &eta=0.9)
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
            //std::cout << "iteration: " << k << ", maxiter: " << maxiter << std::endl;
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
            alpha = compute_decay(dec,k,xk,alpha0,coeff,f,df);

            residual = std::abs( f(xk) - f(xtmp) );
            vec_norm = norm(grad);
        }

        return xk;
    }


    inline Point adam( const Decay &dec,
                            std::function<double(Point)> f, std::function<Point(Point)> df,
                            const Point &x0, const double &alpha0, const unsigned int &maxiter,
                            const double &tols, const double &tolr, const double &coeff=0.9,
                            const double &beta1=0.9, const double &beta2=0.999,
                            const double &epsilon = 1.0e-8)

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
        * @param beta1    rate for first moment
        * @param beta2    rate for second moment
        * @param epsilon  regularizer for denominator 
        * @return res     double containing the optimization point
        **/
    {
        Point m; // first moment
        Point v; // second moment

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
            //std::cout << "iteration: " << k << ", maxiter: " << maxiter << std::endl;
            // update gradient
            grad = df(xk);

            // momentum method element-wise
            for (size_t i=0; i<xk.size(); ++i)
            {
                // update momentum
                m[i] = beta1*m[i] - (1-beta1)*grad[i];
                v[i] = beta2*v[i] - (1-beta2)*grad[i]*grad[i];

                //correct for bias towards 0
                m[i] = m[i]/(1-beta1);
                v[i] = v[i]/(1-beta2);

                // save old point
                xtmp[i] = xk[i];
                // find new point
                //! sqrt(v) length of vector or element-wise?
                xk[i] -= alpha * m[i]/(std::sqrt(v[i] + epsilon));
            }
            //print(xk);

            // learning rate
            alpha = compute_decay(dec,k,xk,alpha0,coeff,f,df);

            residual = std::abs( f(xk) - f(xtmp) );
            vec_norm = norm(grad);
        }

        return xk;
    }

    // template <class Function, class Dfunction>
    // Point optimize(Function const &f, Dfunction const &df,
    inline Point optimize(Parameters const &param, Method m, Decay d)
        /**
        * Minimize using a chosen method among ... with adaptive learning rate
        * 
        * @param param    struct with all parameters
        * @param m        optimization method
        * @param d        adaptive lerning rate method
        * 
        * @return res     double containing the optimization point
        **/
    {
        Point res(param.x0.size());
        
        std::cout << "Optimization strategy: ";
        switch(m)
        {
            case gd : 
                res = gradient_descent(d, param.f,param.df,param.x0,
                param.alpha0,param.maxiter,param.tols,param.tolr,
                param.coeff);
                std::cout << "Gradient descent method" << std::endl;
                break;
            case moment : 
                res = momentum(d, param.f,param.df,param.x0,param.alpha0,
                param.maxiter,param.tols,param.tolr, param.coeff,
                param.eta);
                std::cout << "Momentum method" << std::endl;
                break;
            case ad :
                res = adam(d, param.f,param.df,param.x0,param.alpha0,
                param.maxiter,param.tols,param.tolr, param.coeff,
                param.eta, param.beta2, param.eps);
                std::cout << "Momentum method" << std::endl;
                break;
        }

        std::cout << "Adaptive learning rate method: ";
        switch (d) {
            case arm :
                std::cout << "Armijo" << std::endl;
                break;
            case inv :
                std::cout << "Inverse decay" << std::endl;
                break;
            case exp :
                std::cout << "Exponential decay" << std::endl;
                break;
        
        }
        return res;
    }


} // namespace optimization

#endif