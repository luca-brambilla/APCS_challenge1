// Luca Brambilla
// 10510718 - 919812


#include <muParser.h>
#include <memory>
#include <string>

class MuparserFun2D
{
public:
  MuparserFun2D(const MuparserFun2D &m)
    : m_parser(m.m_parser)
  {
    m_parser.DefineVar("t", &m_varx);
    m_parser.DefineVar("y", &m_varz);

  };

  MuparserFun2D(const std::string &s)
  {
    try
      {
        m_parser.DefineVar("t", &m_varx);
        m_parser.DefineVar("y", &m_varz);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &z)
  {
    double y = 0;

    m_varx = x;
    m_varz = z;

    try
      {
        y = m_parser.Eval();
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
    return y;
  };

private:
  double     m_varx;
  double     m_varz;
  mu::Parser m_parser;
};



class MuparserFun1D
{
public:
  MuparserFun1D(const MuparserFun1D &m)
    : m_parser(m.m_parser)
  {
    m_parser.DefineVar("t", &m_var);
  };

  MuparserFun1D(const std::string &s)
  {
    try
      {
        m_parser.DefineVar("t", &m_var);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x)
  {
    double y = 0;

    m_var = x;

    try
      {
        y = m_parser.Eval();
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
    return y;
  };

private:
  double     m_var;
  mu::Parser m_parser;
};