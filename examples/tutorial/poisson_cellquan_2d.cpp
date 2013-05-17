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


/** @brief A tag class used for storing the permittivity (i.e. a cell quantity) with ViennaData */
struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & /*other*/) const { return false; }
};

int main()
{
  typedef viennagrid::config::triangular_2d                             ConfigType;
  typedef viennagrid::result_of::domain<ConfigType>::type               DomainType;
  typedef viennagrid::segment_t<ConfigType>                             SegmentType;

  typedef viennagrid::result_of::ncell_range<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  typedef viennagrid::result_of::ncell_range<SegmentType, 2>::type   CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type           CellIterator;
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
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(my_domain, "../examples/data/square224.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    return EXIT_FAILURE;
  }
  
  
  //
  // Specify Poisson equation with inhomogeneous permittivity:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol v(0, viennamath::test_tag<>());   //an unknown function used for PDE specification
  viennafem::cell_quan<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( permittivity_key() );  

  //the strong form (not yet functional because of ViennaMath limitations)
  //Equation poisson_equ = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), 0);

  //the weak form:
  Equation poisson_equ = viennamath::make_equation( 
                          viennamath::integral(viennamath::symbolic_interval(),
                                               permittivity * (viennamath::grad(u) * viennamath::grad(v)) ),
                          0);

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
    // Boundary condition: 0 at left boundary, 1 at right boundary
    if ( (*vit)[0] == 0.0)
      viennafem::set_dirichlet_boundary(*vit, 0.0);
    else if ( (*vit)[0] == 1.0)
      viennafem::set_dirichlet_boundary(*vit, 1.0);
    
  }
  
  
  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_assembler fem_assembler;

  
  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  for (size_t i=0; i<my_domain.segments().size(); ++i)
  {
    //set permittivity:
    CellContainer cells = viennagrid::ncells<2>(my_domain.segments()[i]);
    for (CellIterator cit  = cells.begin();
                      cit != cells.end();
                    ++cit)
    {
      if (i==0) //Si
        viennadata::access<permittivity_key, double>(permittivity_key())(*cit) = 3.9; 
      else //SiO2
        viennadata::access<permittivity_key, double>(permittivity_key())(*cit) = 11.9; 
    }
    
    
    fem_assembler(viennafem::make_linear_pde_system(poisson_equ, 
                                                    u,
                                                    viennafem::make_linear_pde_options(0, 
                                                                                       viennafem::lagrange_tag<1>(),
                                                                                       viennafem::lagrange_tag<1>())
                                                  ),
                  my_domain.segments()[i],
                  system_matrix,
                  load_vector
                );
  }
  
  VectorType pde_result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::cg_tag());

  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result, "poisson_cellquan_2d", my_domain, 0);
  
  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
