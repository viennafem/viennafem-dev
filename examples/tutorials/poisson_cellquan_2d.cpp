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


/** @brief A tag class used for storing the permittivity (i.e. a cell quantity) with ViennaData */
struct permittivity_key
{
  // Operator< is required for compatibility with std::map
  bool operator<(permittivity_key const & /*other*/) const { return false; }
};

int main()
{
  typedef viennagrid::triangular_2d_mesh                                                  DomainType;
  typedef viennagrid::result_of::segmentation<DomainType>::type                           SegmentationType;
  typedef SegmentationType::iterator                                                      SegmentationIteratorType;
  typedef viennagrid::result_of::segment_handle<SegmentationType>::type                   SegmentType;
  typedef viennagrid::result_of::cell_tag<DomainType>::type                               CellTagType;
  typedef viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type        VertexType;
  typedef viennagrid::result_of::element<DomainType, CellTagType>::type                   CellType;

  typedef viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type  VertexContainerType;
  typedef viennagrid::result_of::iterator<VertexContainerType>::type                      VertexIteratorType;
  typedef viennagrid::result_of::element_range<SegmentType, CellTagType>::type            CellOnSegmentContainerType;
  typedef viennagrid::result_of::iterator<CellOnSegmentContainerType>::type               CellOnSegmentIteratorType;

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
    viennagrid::io::netgen_reader my_netgen_reader;
    my_netgen_reader(my_domain, segments, "../examples/data/square224.mesh");
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
  viennafem::cell_quan<CellType, viennamath::expr::interface_type>  permittivity; permittivity.wrap_constant( storage, permittivity_key() );

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
  VertexContainerType vertices = viennagrid::elements<VertexType>(my_domain);
  for (VertexIteratorType vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    // Boundary condition: 0 at left boundary, 1 at right boundary
    if ( viennagrid::point(my_domain, *vit)[0] == 0.0)
      viennafem::set_dirichlet_boundary(storage, *vit, 0.0);
    else if ( viennagrid::point(my_domain, *vit)[0] == 1.0)
      viennafem::set_dirichlet_boundary(storage, *vit, 1.0);

  }


  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_assembler<StorageType> fem_assembler(storage);


  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  std::size_t si = 0;
  for(SegmentationIteratorType sit = segments.begin(); sit != segments.end(); sit++)
  {
    //set permittivity:
    CellOnSegmentContainerType cells = viennagrid::elements<CellType>(*sit);
    for (CellOnSegmentIteratorType cit  = cells.begin();
                                   cit != cells.end();
                                 ++cit)
    {
      if (si == 0) //Si
        viennadata::access<permittivity_key, double>(storage, permittivity_key(), *cit) = 3.9;
      else //SiO2
        viennadata::access<permittivity_key, double>(storage, permittivity_key(), *cit) = 11.9;
    }


    fem_assembler(viennafem::make_linear_pde_system(poisson_equ,
                                                    u,
                                                    viennafem::make_linear_pde_options(0,
                                                                                       viennafem::lagrange_tag<1>(),
                                                                                       viennafem::lagrange_tag<1>())
                                                  ),
                  *sit,
                  system_matrix,
                  load_vector
                );
  }

  VectorType pde_result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::cg_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, pde_result) - load_vector) << std::endl;

  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result, "poisson_cellquan_2d", my_domain, segments, storage, 0);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
