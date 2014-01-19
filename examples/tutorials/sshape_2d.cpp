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
  typedef viennagrid::triangular_2d_mesh                                                  DomainType;
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
    my_reader(my_domain, segments, "../examples/data/sshape2d.mesh");
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
  VertexContainer vertices = viennagrid::elements<VertexType>(my_domain);
  for (VertexIterator vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    //boundary for first equation: Homogeneous Dirichlet everywhere
    if (viennagrid::point(my_domain, *vit)[1] == 3.0 || viennagrid::point(my_domain, *vit)[0] == 5.0 )
      viennafem::set_dirichlet_boundary(storage, *vit, 0.0);

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
                system_matrix,
                load_vector
               );

  VectorType pde_result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::cg_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, pde_result) - load_vector) << std::endl;

  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result, "sshape_2d", my_domain, segments, storage, 0);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
