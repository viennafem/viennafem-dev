#ifndef VIENNAFEM_LOG_LATEX_HPP
#define VIENNAFEM_LOG_LATEX_HPP

/* =========================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */


#include <vector>
#include <fstream>
#include "viennafem/forwards.h"
#include "viennafem/log/interface.hpp"
#include "viennafem/cell_quan.hpp"

#include "viennamath/manipulation/latex.hpp"

/** @file  viennafem/log/latex.hpp
    @brief Defines a logger writing a LaTeX document
*/

namespace viennafem
{
  template <typename InterfaceType>
  class latex_logger;

  namespace detail
  {

    template <typename EquationArray, typename InterfaceType>
    void write_strong_form(EquationArray const & pde_system,
                          latex_logger<InterfaceType> & log);

    template <typename EquationArray, typename InterfaceType>
    void write_weak_form(EquationArray const & weak_form,
                        latex_logger<InterfaceType> & log);

    template <typename EquationArray, typename InterfaceType>
    void write_coordinated_weak_form(EquationArray const & weak_form,
                                    latex_logger<InterfaceType> & log);

    template <typename EquationArray, typename InterfaceType>
    void write_transformed_weak_form(EquationArray const & weak_form,
                                    latex_logger<InterfaceType> & log);

    template <typename EquationArray, typename InterfaceType>
    void write_test_and_trial_space(EquationArray const & test_space,
                                    EquationArray const & trial_space,
                                    latex_logger<InterfaceType> & log);

    template <typename InterfaceType>
    void write_linear_solver_stats(latex_logger<InterfaceType> & log);

  }


  //
  // cell_quan_handler
  //
  /** @brief Defines a custom LaTeX processor for cell_quan expressions
   *
   * @tparam CellType
   * @tparam InterfaceType
   */
  template <typename CellType, typename StorageType, typename InterfaceType>
  class rt_latex_dt_dx_processor : public viennamath::rt_latex_processor_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type              NumericType;
      typedef viennafem::cell_quan<CellType, InterfaceType>     CellQuanType;

    public:
      viennamath::rt_latex_processor_interface<InterfaceType> * clone() const { return new rt_latex_dt_dx_processor(); }

      std::string process(InterfaceType const * ptr, bool use_parenthesis, viennamath::rt_latex_translator<InterfaceType> const & translator) const
      {
        if (dynamic_cast< const CellQuanType * >(ptr) != NULL)
        {
          const CellQuanType * temp = dynamic_cast< const CellQuanType * >(ptr);
          return process_impl(*temp, use_parenthesis, translator);
        }

        return "";
      }

      bool customize(InterfaceType const * /*e*/, std::string const & /*str*/)
      {
        return false;
      }

    private:

      std::string process_impl(CellQuanType const & e, bool /*use_parenthesis*/, viennamath::rt_latex_translator<InterfaceType> const & /*translator*/) const
      {
        typedef viennamath::expr      ExpressionType;

        std::stringstream ss;

        viennamath::variable x(0);
        viennamath::variable y(1);
        viennamath::variable z(2); 

        std::auto_ptr< detail::cell_quan_interface<CellType> > temp(e.wrapper().clone());  //create a clone in order to dispatch with respect to the type

        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,0>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,0>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\xi}{\\partial x} ";
        }
        else 
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,0>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,0>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\eta}{\\partial x} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,0>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,0>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\nu}{\\partial x} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,1>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,1>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\xi}{\\partial y} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,1>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,1>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\eta}{\\partial y} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,1>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,1>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\nu}{\\partial y} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,2>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<0,2>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\xi}{\\partial z} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,2>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<1,2>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\eta}{\\partial z} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,2>, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::dt_dx_key<2,2>, NumericType, CellType>::type > * >(temp.get()) != NULL)
        {
          ss << " \\frac{\\partial \\nu}{\\partial z} ";
        }
        else
        if (dynamic_cast< detail::cell_quan_expr<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::det_dF_dt_key, ExpressionType, CellType>::type > * >(temp.get()) != NULL
            || dynamic_cast< detail::cell_quan_constant<CellType, typename viennadata::result_of::accessor<StorageType, viennafem::det_dF_dt_key, NumericType, CellType>::type > * >(temp.get()) != NULL)        {
          ss << " \\mathrm{det} F ";
        }
        else
          std::cerr << "Warning: could not find string for cell_quan!" << std::endl;





//        if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<0,0>, ExpressionType> * >(temp.get()) != NULL
//            || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<0,0>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\xi}{\\partial x} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<1,0>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<1,0>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\eta}{\\partial x} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<2,0>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<2,0>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\nu}{\\partial x} ";
//        }

//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<0,1>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<0,1>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\xi}{\\partial y} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<1,1>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<1,1>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\eta}{\\partial y} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<2,1>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<2,1>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\nu}{\\partial y} ";
//        }

//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<0,2>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<0,2>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\xi}{\\partial z} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<1,2>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<1,2>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\eta}{\\partial z} ";
//        }
//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::dt_dx_key<2,2>, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::dt_dx_key<2,2>, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\frac{\\partial \\nu}{\\partial z} ";
//        }

//        else if (dynamic_cast< detail::cell_quan_expr<CellType, viennafem::det_dF_dt_key, ExpressionType> * >(temp.get()) != NULL
//                 || dynamic_cast< detail::cell_quan_constant<CellType, viennafem::det_dF_dt_key, NumericType> * >(temp.get()) != NULL)
//        {
//          ss << " \\mathrm{det} F ";
//        }
//        else
//          std::cerr << "Warning: could not find string for cell_quan!" << std::endl;

        return ss.str();
      }

  };



  /** @brief The LaTeX logger class.
   *
   * @tparam InterfaceType    The ViennaMath runtime interface
   */
  template <typename InterfaceType>
  class latex_logger : public logger_interface<InterfaceType>
  {
      typedef logger_interface<InterfaceType>    BaseType;
      typedef typename BaseType::EquationType             EquationType;

    public:
      latex_logger(std::string const & filename) : stream_(filename.c_str()), filename_(filename), something_written(false)
      {
        viennamath::function_symbol u0(0, viennamath::unknown_tag<>());
        viennamath::function_symbol u1(1, viennamath::unknown_tag<>());
        viennamath::function_symbol u2(2, viennamath::unknown_tag<>());

        viennamath::function_symbol v0(0, viennamath::test_tag<>());
        viennamath::function_symbol v1(1, viennamath::test_tag<>());
        viennamath::function_symbol v2(2, viennamath::test_tag<>());

        translator_.customize(u0, "u_0");
        translator_.customize(u1, "u_1");
        translator_.customize(u2, "u_2");

        translator_.customize(v0, "v_0");
        translator_.customize(v1, "v_1");
        translator_.customize(v2, "v_2");
      }

      latex_logger(const latex_logger & other) : stream_(other.filename_.c_str()),
                                                 filename_(other.filename_),
                                                 something_written(false),
                                                 translator_(other.translator_) {}

      ~latex_logger()
      {
        //write footer:
        if (something_written)
          stream_ << "\\end{document}" << std::endl;

        stream_.close();
      }

      void write_strong_form(std::vector<EquationType> const & pdes) { viennafem::detail::write_strong_form(pdes, *this); }
      void write_weak_form(std::vector<EquationType> const & pdes) { viennafem::detail::write_weak_form(pdes, *this); }
      void write_coordinated_weak_form(std::vector<EquationType> const & pdes) { viennafem::detail::write_coordinated_weak_form(pdes, *this); }
      void write_transformed_weak_form(std::vector<EquationType> const & pdes) { viennafem::detail::write_transformed_weak_form(pdes, *this); }
      void write_test_and_trial_space(std::vector<viennamath::expr> const & test_space,
                                      std::vector<viennamath::expr> const & trial_space)
      {
        viennafem::detail::write_test_and_trial_space(test_space, trial_space, *this);
      }


      void write_header()
      {
        stream_ << "\\documentclass[11pt]{article}\n";
        stream_ << "\\usepackage{amsmath}\n";
        stream_ << "\\usepackage[cm]{fullpage}\n";
        stream_ << "\\title{ViennaFEM Protocol}\n";
        stream_ << "\n";
        stream_ << "\\begin{document}\n";
        stream_ << "\\begin{center}\n";
        stream_ << "\\LARGE ViennaFEM Protocol\n";
        stream_ << "\\end{center}\n";
        stream_ << "\n";

      }

      template <typename T>
      void stream(T const & t)
      {
        if (!something_written)
        {
          write_header();
          something_written = true;
        }
        stream_ << t;
      }
      void flush() {  }

      viennamath::latex_translator & translator() { return translator_; }

    private:
      latex_logger & operator=(latex_logger const & other);

      std::ofstream stream_;
      const std::string filename_;
      bool something_written;
      viennamath::latex_translator translator_;
  };

  /** @brief Convenience overload for streaming text to a LaTeX logger */
  template <typename InterfaceType, typename T>
  latex_logger<InterfaceType> & operator<<(latex_logger<InterfaceType> & logger, T const & t)
  {
     logger.stream(t);
     return logger;
  }



  namespace detail
  {

    /** @brief Implementation for writing the strong form to LaTeX */
    template <typename EquationArray, typename InterfaceType>
    void write_strong_form(EquationArray const & pdes,
                          latex_logger<InterfaceType> & log)
    {
      log << "\\section{Strong Formulation}\n";
      log << "The strong formulation of the problem is to solve\n";
      log << "\\begin{align}\n";
      log << log.translator()(pdes[0]) << " \n";
      log << "\\end{align}\n";
      log << "in $\\Omega$ with suitable boundary conditions.\n";
    }

    /** @brief Implementation for writing the weak form to LaTeX */
    template <typename EquationArray, typename InterfaceType>
    void write_weak_form(EquationArray const & weak_form,
                        latex_logger<InterfaceType> & log)
    {
      log << "\\section{Weak Formulation}\n";
      log << "After multiplication with a test function $v$, integration over $\\Omega$ and partial integration, the weak formulation of the problem is obtained as\n";
      log << "\\begin{align}\n";
      for (typename EquationArray::const_iterator it = weak_form.begin();
                                                  it != weak_form.end();
                                                ++it)
        log << log.translator()(*it) << " \n";
      log << "\\end{align}\n";
      log << "for all test functions $v$.\n";
    }

    /** @brief Implementation for writing the coordinated weak form to LaTeX */
    template <typename EquationArray, typename InterfaceType>
    void write_coordinated_weak_form(EquationArray const & weak_form,
                                    latex_logger<InterfaceType> & log)
    {
      log << "Written in coordinates for the underlying simulation domain, the weak form is \n";
      log << "\\begin{align}\n";
      for (typename EquationArray::const_iterator it = weak_form.begin();
                                                  it != weak_form.end();
                                                ++it)
        log << log.translator()(*it) << " \n";
      log << "\\end{align}\n";
      log << "for all testfunctions $v$.\n";
    }

    /** @brief Implementation for writing the transformed weak form to LaTeX */
    template <typename EquationArray, typename InterfaceType>
    void write_transformed_weak_form(EquationArray const & weak_form,
                                    latex_logger<InterfaceType> & log)
    {
      // Restriction to the discrete space
      log << "\\section{Discrete Approximation by Finite Elements}\n";
      log << "The discrete approximation to the continuous problem is obtained by a Galerkin approach:\n";
      log << "\\begin{align}\n";
      log << " u_h = \\sum_j \\alpha_j \\varphi_j  \n";
      log << "\\end{align}\n";
      log << "with trial functions $\\varphi_j$.\n";
      //log << "Thus, instead of solving for the continuous function $u$, only the coefficients $\\alpha_i$ need to be computed.\n";
      log << "In a Galerkin approach, test functions $v$ are also chosen from a finite-dimensional space:\n";
      log << "\\begin{align}\n";
      log << " v_h = \\sum_i \\beta_i \\psi_i  \n";
      log << "\\end{align}\n";
      log << "Due to integral transformations, it is sufficient to define the trial and test functions on the reference cell.\n";
      log << "After transformation to the reference cell, the weak form on a cell reads\n";

      // give new name to local variables:
      viennamath::variable xi(0);
      viennamath::variable eta(1);
      viennamath::variable nu(2);

      log.translator().customize(xi, "\\xi");
      log.translator().customize(eta, "\\eta");
      log.translator().customize(nu, "\\nu");

      viennamath::function_symbol u0(0, viennamath::unknown_tag<>());
      viennamath::function_symbol u1(1, viennamath::unknown_tag<>());
      viennamath::function_symbol u2(2, viennamath::unknown_tag<>());

      viennamath::function_symbol v0(0, viennamath::test_tag<>());
      viennamath::function_symbol v1(1, viennamath::test_tag<>());
      viennamath::function_symbol v2(2, viennamath::test_tag<>());

      log.translator().customize(u0, "\\tilde{u}_0");
      log.translator().customize(u1, "\\tilde{u}_1");
      log.translator().customize(u2, "\\tilde{u}_2");

      log.translator().customize(v0, "\\tilde{v}_0");
      log.translator().customize(v1, "\\tilde{v}_1");
      log.translator().customize(v2, "\\tilde{v}_2");

      // [JW] switched to $ math expressions, otherwise the equation is cut
      // off at the right page border ...
      //
      //log << "\\begin{align}\n";
      log << "\\newline\\newline$\n";
      for (typename EquationArray::const_iterator it = weak_form.begin();
                                                  it != weak_form.end();
                                                ++it)
        log << log.translator()(*it) << " \\  . \n";
      //log << "\\end{align}\n";
      log << "$\\newline\\newline\n";


      //
      //log << "write transformed weak form\n";
    }

    /** @brief Implementation for writing the test and trial spaces to LaTeX */
    template <typename EquationArray, typename InterfaceType>
    void write_test_and_trial_space(EquationArray const & test_space,
                                    EquationArray const & trial_space,
                                    latex_logger<InterfaceType> & log)
    {
      log << "The trial functions used for the discrete approximation of $u_h$ on the reference element are \n";
      log << "\\begin{align*}\n";
      std::size_t func_index = 0;
      for (typename EquationArray::const_iterator it = trial_space.begin();
                                                  it != trial_space.end();
                                                ++it)
      {
        if (func_index > 0) //finish previous line
          log << " \\ , \\\\ \n";

        log << "\\varphi_{" << func_index << "} &= " << log.translator()(*it);
        ++func_index;
      }
      log << " \\ . \n \\end{align*}\n";
      log << "The test functions on the reference element are\n";
      log << "\\begin{align*}\n";
      func_index = 0;
      for (typename EquationArray::const_iterator it = test_space.begin();
                                                  it != test_space.end();
                                                ++it)
      {
        if (func_index > 0) //finish previous line
          log << " \\ , \\\\ \n";

        log << "\\psi_{" << func_index << "} &= " << log.translator()(*it);
        ++func_index;
      }
      log << " \\ . \n \\end{align*}\n";
      log << "\n";
    }

    /** @brief Implementation for writing linear solver statistics to LaTeX */
    template <typename InterfaceType>
    void write_linear_solver_stats(latex_logger<InterfaceType> & log)
    {
      log << "write linear solver stats\n";
    }

  } //namespace detail

} //namespace viennafem

#endif
