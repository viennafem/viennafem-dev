/* =========================================================================
   Copyright (c) 2012-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */

// remove assert() statements and the like in order to get reasonable performance
#ifndef NDEBUG
  #define NDEBUG
#endif

// include necessary system headers
#include <iostream>

// ViennaFEM includes:
#include "viennafem/fem.hpp"
#include "viennafem/io/vtk_writer.hpp"

// ViennaGrid includes:
#include "viennagrid/forwards.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"

// Boost.uBLAS includes:
#include <boost/numeric/ublas/io.hpp>
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
  typedef viennagrid::hexahedral_3d_mesh                                                  DomainType;
  typedef viennagrid::result_of::segmentation<DomainType>::type                           SegmentationType;

  typedef viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type        VertexType;
  typedef viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type  VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type                          VertexIterator;

  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;

  //
  // Create a domain from file
  //
  DomainType my_domain;
  SegmentationType segments(my_domain);

  //
  // Create a storage object
  //
  typedef viennadata::storage<> StorageType;
  StorageType   storage;

  try
  {
    viennagrid::io::netgen_reader my_reader;
    my_reader(my_domain, segments, "../examples/data/cube343_hex.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    exit(EXIT_FAILURE);
  }


  //
  // Specify two PDEs:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  Equation poisson_equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  Equation poisson_equ_2 = viennamath::make_equation( viennamath::laplace(u), -1);

  MatrixType system_matrix_1, system_matrix_2;
  VectorType load_vector_1, load_vector_2;

  //
  // Setting boundary information on domain (this should come from device specification)
  //
  //setting some boundary flags:
  VertexContainer vertices = viennagrid::elements<VertexType>(my_domain);
  for (VertexIterator vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    // First equation: Homogeneous boundary conditions at x=0, x=1, y=0, or y=1
    if ( viennagrid::point(my_domain, *vit)[0] == 0.0 || viennagrid::point(my_domain, *vit)[0] == 1.0
         || viennagrid::point(my_domain, *vit)[1] == 0.0 || viennagrid::point(my_domain, *vit)[1] == 1.0 )
      viennafem::set_dirichlet_boundary(storage, *vit, 0.0, 0);  //simulation with ID 0 uses homogeneous boundary data

    // Boundary for second equation (ID 1): 0 at left boundary, 1 at right boundary
    if ( viennagrid::point(my_domain, *vit)[0] == 0.0)
      viennafem::set_dirichlet_boundary(storage, *vit, 0.0, 1);
    else if ( viennagrid::point(my_domain, *vit)[0] == 1.0)
      viennafem::set_dirichlet_boundary(storage, *vit, 1.0, 1);
  }


  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_assembler<StorageType> fem_assembler(storage);


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
                system_matrix_1,
                load_vector_1
               );

  fem_assembler(viennafem::make_linear_pde_system(poisson_equ_2,
                                                  u,
                                                  viennafem::make_linear_pde_options(1,
                                                                                     viennafem::lagrange_tag<1>(),
                                                                                     viennafem::lagrange_tag<1>())
                                                 ),
                my_domain,
                system_matrix_2,
                load_vector_2
               );

  VectorType pde_result_1 = viennacl::linalg::solve(system_matrix_1, load_vector_1, viennacl::linalg::cg_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix_1, pde_result_1) - load_vector_1) << std::endl;

  VectorType pde_result_2 = viennacl::linalg::solve(system_matrix_2, load_vector_2, viennacl::linalg::cg_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix_2, pde_result_2) - load_vector_2) << std::endl;


  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result_1, "poisson_3d_hex_1", my_domain, segments, storage, 0);
  viennafem::io::write_solution_to_VTK_file(pde_result_2, "poisson_3d_hex_2", my_domain, segments, storage, 1);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
