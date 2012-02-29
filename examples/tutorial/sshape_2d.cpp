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

// Boost.uBLAS includes:
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>


//ViennaCL includes:
#ifndef VIENNACL_HAVE_UBLAS
 #define VIENNACL_HAVE_UBLAS
#endif
    
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"


int main()
{
  typedef viennagrid::config::triangular_2d                             ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type         DomainType;

  typedef viennagrid::result_of::ncell_range<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  typedef viennagrid::result_of::ncell<ConfigType, 2>::type              CellType;
  
  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;
  
  typedef viennafem::boundary_key      BoundaryKey;
  
  //
  // Create a domain from file
  //
  DomainType my_domain;
  
  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(my_domain, "../examples/data/sshape2d.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }
  
  
  //
  // Specify PDE:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  
  MatrixType system_matrix;
  VectorType load_vector;
  
  //
  // Setting boundary information on domain (this should come from device specification)
  //
  //setting some boundary flags:
  VertexContainer vertices = viennagrid::ncells<0>(my_domain);
  for (VertexIterator vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    //boundary for first equation: Homogeneous Dirichlet everywhere
    if (vit->point()[1] == 3.0 || vit->point()[0] == 5.0 )
      viennafem::set_dirichlet_boundary(*vit, 0.0);
    
  }
  
  
  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_assembler fem_assembler;

  
  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  fem_assembler(viennafem::make_linear_pde_system(poisson_equ_1, 
                                                  u,
                                                  viennafem::make_linear_pde_options(0, 
                                                                                     viennafem::lagrange_tag<1>(),
                                                                                     viennafem::lagrange_tag<1>())
                                                 ),
                my_domain,
                system_matrix,
                load_vector
               );
  
  VectorType pde_result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::cg_tag());
  
  
  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result, "sshape_2d", my_domain, 0);
  
  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
