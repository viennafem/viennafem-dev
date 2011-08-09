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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>


//ViennaCL includes:
#ifndef VIENNACL_HAVE_UBLAS
 #define VIENNACL_HAVE_UBLAS
#endif
    
#ifdef USE_OPENCL
  #include "viennacl/matrix.hpp"
  #include "viennacl/vector.hpp"
#endif
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"


//      
// Solve system of linear equations:
//
template <typename MatrixType, typename VectorType>
VectorType solve(MatrixType const & system_matrix,
                 VectorType const & load_vector)
{
  typedef typename VectorType::value_type        numeric_type;
  VectorType result(load_vector.size());
  
  std::cout << "* solve(): Solving linear system" << std::endl;

#ifdef USE_OPENCL
  viennacl::matrix<viennafem::numeric_type> vcl_matrix(load_vector.size(), load_vector.size());
  viennacl::vector<viennafem::numeric_type> vcl_rhs(load_vector.size());
  viennacl::vector<viennafem::numeric_type> vcl_result(load_vector.size());
  
  viennacl::copy(system_matrix, vcl_matrix);
  viennacl::copy(load_vector, vcl_rhs);
  
  vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::cg_tag());
  
  viennacl::copy(vcl_result, result);
#else
  result = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::cg_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, result) - load_vector) << std::endl;
#endif
    
  //std::cout << load_vector << std::endl;
  
  //print solution:
  //std::cout << "Solution: ";
  //for (size_t i=0; i<ublas_result.size(); ++i)
  //  std::cout << ublas_result(i) << " ";
  //std::cout << std::endl;
  //std::cout << std::endl;

  return result;
}

struct permittivity_key
{
  bool operator<(permittivity_key const & other) const { return false; }
};

int main()
{
  typedef viennagrid::config::triangular_2d                             ConfigType;
  typedef viennagrid::domain<ConfigType>                                DomainType;
  typedef viennagrid::segment_t<ConfigType>                             SegmentType;

  typedef viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  typedef viennagrid::result_of::ncell_container<SegmentType, 2>::type   CellContainer;
  typedef viennagrid::result_of::iterator<CellContainer>::type           CellIterator;
  typedef viennagrid::result_of::ncell_type<ConfigType, 2>::type              CellType;
  
  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;

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
  // Specify Poisson equation with inhomogeneous permittivity:
  //
  FunctionSymbol u(0, viennamath::unknown_tag<>());   //an unknown function used for PDE specification
  FunctionSymbol v(0, viennamath::test_tag<>());   //an unknown function used for PDE specification
  viennafem::cell_quan<CellType, viennamath::expr<>::interface_type>  permittivity; permittivity.wrap( permittivity_key() );  

  //the strong form (not yet functional because of ViennaMath limitations)
  //Equation poisson_equ = viennamath::make_equation( viennamath::div(permittivity * viennamath::grad(u)), 0);

  //the weak form:
  Equation poisson_equ = viennamath::make_equation( 
                          viennamath::integral(viennamath::Omega(),
                                               permittivity * (viennamath::grad(u) * viennamath::grad(v)),
                                               viennamath::symbolic_tag()),
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
    //boundary for second equation: 0 at left boundary, 1 at right boundary
    if (vit->getPoint()[0] == 0.0)
    {
      viennadata::access<BoundaryKey, bool>(BoundaryKey(0))(*vit) = true;
      viennadata::access<BoundaryKey, double>(BoundaryKey(0))(*vit) = 0.0;
    }
    else if (vit->getPoint()[0] == 1.0)
    {
      viennadata::access<BoundaryKey, bool>(BoundaryKey(0))(*vit) = true;
      viennadata::access<BoundaryKey, double>(BoundaryKey(0))(*vit) = 1.0;
    }
    else
      viennadata::access<BoundaryKey, bool>(BoundaryKey(0))(*vit) = false;
    
  }
  
  
  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_solver fem_assembler;

  
  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  for (size_t i=0; i<my_domain.segment_size(); ++i)
  {
    //set permittivity:
    CellContainer cells = viennagrid::ncells<2>(my_domain.segment(i));
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
                                                                                      viennafem::LinearBasisfunctionTag(),
                                                                                      viennafem::LinearBasisfunctionTag())
                                                  ),
                  my_domain.segment(i),
                  system_matrix,
                  load_vector
                );
  }
  
  //std::cout << poisson_config_1.load_vector() << std::endl;
  
  VectorType pde_result = solve(system_matrix, load_vector);

  //std::cout << "RESULT" << std::endl;
  //std::cout << pde_result_1 << std::endl;
  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  viennafem::io::write_solution_to_VTK_file(pde_result, "poisson_1", my_domain, 0);
  
  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
