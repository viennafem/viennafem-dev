#ifndef VIENNAFEM_LOG_LATEX_HPP
#define VIENNAFEM_LOG_LATEX_HPP

/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */


#include <vector>
#include <fstream>
#include "viennafem/forwards.h"
#include "viennafem/log/interface.hpp"

#include "viennamath/manipulation/latex.hpp"

namespace viennafem
{
  template <typename InterfaceType>
  class latex_logger;
  
  
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
      
      void write_strong_form(std::vector<EquationType> const & pdes) { viennafem::write_strong_form(pdes, *this); }
      void write_weak_form(std::vector<EquationType> const & pdes) { viennafem::write_weak_form(pdes, *this); }
      void write_coordinated_weak_form(std::vector<EquationType> const & pdes) { viennafem::write_coordinated_weak_form(pdes, *this); }
      void write_transformed_weak_form(std::vector<EquationType> const & pdes) { viennafem::write_transformed_weak_form(pdes, *this); }
      void write_test_and_trial_space(std::vector<viennamath::expr> const & test_space,
                                      std::vector<viennamath::expr> const & trial_space)
      {
        viennafem::write_test_and_trial_space(test_space, trial_space, *this);
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
  

  template <typename InterfaceType, typename T>
  latex_logger<InterfaceType> & operator<<(latex_logger<InterfaceType> & logger, T const & t)
  {
     logger.stream(t);
     return logger;
  }

  template <typename EquationArray, typename InterfaceType>
  void write_strong_form(EquationArray const & pdes,
                         latex_logger<InterfaceType> & log)
  {
    log << "\\section{Strong Formulation}\n";
    log << "The strong formulation of the problem is to solve\n";
    log << "\\begin{align}\n";
    log << log.translator()(pdes[0]) << " \n";
    log << "\\end{align}\n";
    log << "in $\\Omega$ and equipped with suitable boundary conditions.\n";
  }

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

  template <typename EquationArray, typename InterfaceType>
  void write_transformed_weak_form(EquationArray const & weak_form,
                                   latex_logger<InterfaceType> & log)
  {
    // Restriction to the discrete space
    log << "\\section{Discrete Approximation by Finite Elements}\n";
    log << "The discrete approximation to the continuous problem is obtained by a Galerkin approach:\n";
    log << "\\begin{align}\n";
    log << " u_h = \\sum \\alpha_j \\varphi_j  \n";
    log << "\\end{align}\n";
    log << "with trial functions $\\varphi_j$\n";
    //log << "Thus, instead of solving for the continuous function $u$, only the coefficients $\\alpha_i$ need to be computed.\n";
    log << "In a Galerkin approach, test functions $v$ are also chosen from a finite-dimensional space:\n";
    log << "\\begin{align}\n";
    log << " v_h = \\sum \\beta_i \\psi_i  \n";
    log << "\\end{align}\n";
    log << "Due to integral transformations, it is sufficient to define the trial and test functions on the reference cell.\n";
    log << "After transformation to the reference cell, the weak form on a cell reads\n";

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
    
    log << "\\begin{align}\n";
    for (typename EquationArray::const_iterator it = weak_form.begin();
                                                it != weak_form.end();
                                              ++it)
      log << log.translator()(*it) << " \n";
    log << "\\end{align}\n";
   
    
    //
    //log << "write transformed weak form\n";
  }

  template <typename EquationArray, typename InterfaceType>
  void write_test_and_trial_space(EquationArray const & test_space,
                                  EquationArray const & trial_space,
                                  latex_logger<InterfaceType> & log)
  {
    log << "write test and trial space\n";
  }
  
  template <typename InterfaceType>
  void write_linear_solver_stats(latex_logger<InterfaceType> & log)
  {
    log << "write linear solver stats\n";
  }



}

#endif
