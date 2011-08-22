/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

// include necessary system headers
#include <iostream>

// ViennaFEM includes:
//#include "viennafem/afftrans.hpp"
//#include "viennafem/dtdx_tetrahedron.h"
#include "viennafem/typelist.h"
#include "viennafem/forwards.h"
//#include "viennafem/assembling.hpp"
//#include "viennafem/mapping.hpp"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/unknown_config.hpp"
#include "viennafem/pde_assembler.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/io/vtk_writer.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include <viennagrid/config/simplex.hpp>
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"


struct testkey 
{
  bool operator<(testkey const & other) const { return false; }
};

int main()
{
  typedef viennagrid::config::triangular_2d                             ConfigType;
  typedef viennagrid::domain<ConfigType>         DomainType;

  typedef viennagrid::result_of::ncell_range<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  typedef viennagrid::result_of::ncell<ConfigType, 2>::type              CellType;
  
  typedef viennamath::function_symbol<>   FunctionSymbol;
  typedef viennamath::equation<>          Equation;
  
  typedef viennafem::boundary_key      BoundaryKey;
  
  //
  // Create a domain from file
  //
  DomainType my_domain;

  my_domain.create_segments(2);
  
  try
  {
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(my_domain, "../examples/data/square224.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }
  
  
  //
  // Specify two PDEs:
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
  
  viennafem::cell_quan<CellType, viennamath::expr<>::interface_type>  permittivity; permittivity.wrap( testkey() );  
  Equation equ_6 = viennamath::make_equation( viennamath::div( permittivity * viennamath::grad(u)), -1);
  std::cout << "Strong form: " << equ_6 << std::endl;
  std::cout << "Weak form: " << viennafem::make_weak_form(equ_6) << std::endl;  
  std::cout << "-------------" << std::endl;
  
  std::cout << "*************************************" << std::endl;
  std::cout << "* Weak forms finished successfully! *" << std::endl;
  std::cout << "*************************************" << std::endl;
  return EXIT_SUCCESS;
}
