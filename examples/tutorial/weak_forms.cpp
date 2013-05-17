/* =======================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaMath - Symbolic and Numerical Math in C++
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaMath base directory
======================================================================= */


// include necessary system headers
#include <iostream>

// ViennaFEM includes:
#include "viennafem/fem.hpp"
#include "viennafem/io/vtk_writer.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include "viennagrid/config/simplex.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"

/*
 *   Tutorial: Demonstration of weak form manipulation inside ViennaFEM
 */


// A test key used later for storing a cell quantity
struct testkey 
{
  bool operator<(testkey const & /*other*/) const { return false; }
};

int main()
{
  typedef viennagrid::config::triangular_2d                       ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type         DomainType;
  typedef viennagrid::result_of::ncell<ConfigType, 2>::type       CellType;
  
  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;
  
  //
  // Specify PDEs and derive weak form:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  std::cout << "Strong form: " << equ_1 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_1) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  Equation equ_2 = viennamath::make_equation( viennamath::laplace(u) + u, -1);
  std::cout << "Strong form: " << equ_2 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_2) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  Equation equ_3 = viennamath::make_equation( viennamath::div( 3.0 * viennamath::grad(u)), -1);
  std::cout << "Strong form: " << equ_3 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_3) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  Equation equ_4 = viennamath::make_equation( viennamath::div( 4.0 * viennamath::grad(u)) + u, -1);
  std::cout << "Strong form: " << equ_4 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_4) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  Equation equ_5 = viennamath::make_equation( u, -1);
  std::cout << "Strong form: " << equ_5 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_5) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  viennafem::cell_quan<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( testkey() );  
  Equation equ_6 = viennamath::make_equation( viennamath::div( permittivity * viennamath::grad(u)), -1);
  std::cout << "Strong form: " << equ_6 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_6) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  std::cout << "*************************************" << std::endl;
  std::cout << "* Weak forms finished successfully! *" << std::endl;
  std::cout << "*************************************" << std::endl;
  return EXIT_SUCCESS;
}
