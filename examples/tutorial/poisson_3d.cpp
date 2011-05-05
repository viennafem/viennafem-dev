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
#include "viennafem/dtdx_triangle.h"
#include "viennafem/typelist.h"
#include "viennafem/forwards.h"
#include "viennafem/BFStock.hpp"
//#include "viennafem/assembling.hpp"
//#include "viennafem/mapping.hpp"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/eval.hpp"
#include "viennafem/unknown_config.hpp"
#include "viennafem/pde_solver.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include <viennagrid/config/simplex.hpp>
#include "viennagrid/io/sgf_reader.hpp"
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
// Configuration class for a triangular domain
//
struct TetrahedronConfig
{
  typedef double                                    numeric_type;
  typedef viennagrid::three_dimensions_tag          dimension_tag;
  typedef viennagrid::tetrahedron_tag               cell_tag;

  //multigrid:
  //typedef viennagrid::full_multigrid_tag                       multigrid_tag;
  typedef viennagrid::no_multigrid_tag             multigrid_tag;
};


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

template <typename VectorType,
          typename DomainType,
          typename PDEConfig>
void write_solution_to_VTK_file(VectorType const & result,
                                std::string filename,
                                DomainType const & domain,
                                PDEConfig const & config)
{
  typedef typename viennagrid::result_of::const_ncell_container<DomainType, 0>::type    VertexContainer;
  typedef typename viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  
  typedef typename PDEConfig::mapping_key_type          MappingType;
  typedef typename PDEConfig::boundary_key_type         BoundaryType;
  
  std::cout << "* write_solution_to_VTK_file(): Writing result on mesh for later export" << std::endl;
  VertexContainer vertices = viennagrid::ncells<0>(domain);
  for (VertexIterator vit = vertices.begin();
       vit != vertices.end();
       ++vit)
  {
     long cur_index = viennadata::access<MappingType, long>(config.mapping_key())(*vit);
     if (cur_index > -1)
       viennadata::access<std::string, double>("vtk_data")(*vit) = result[cur_index];
     else //use Dirichlet boundary data:
       viennadata::access<std::string, double>("vtk_data")(*vit) = 
        viennadata::access<BoundaryType, double>(config.boundary_key())(*vit);
  }

  std::cout << "* write_solution_to_VTK_file(): Writing data to '"
            << filename
            << "' (can be viewed with e.g. Paraview)" << std::endl;

  viennagrid::io::Vtk_writer<DomainType> my_vtk_writer;
  my_vtk_writer.writeDomain(domain, filename);  
}



int main()
{
  typedef viennagrid::config::tetrahedral_3d                             ConfigType;
  typedef viennagrid::domain<ConfigType>         DomainType;

  typedef viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  typedef viennagrid::result_of::ncell_type<ConfigType, 3>::type              CellType;
  
  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;
  
  typedef viennamath::function_symbol<viennafem::fem_expression_interface<viennafem::numeric_type, CellType> >   FunctionSymbol;
  typedef viennamath::equation<viennafem::fem_expression_interface<viennafem::numeric_type, CellType> >          Equation;
  
  //
  // Create a domain from file
  //
  DomainType my_domain;
  
  try
  {
    viennagrid::io::sgf_reader my_sgf_reader;
    my_sgf_reader(my_domain, "../examples/data/cube3072.sgf");
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
  
  //
  // Create PDE config (where to find boundary information, where to store mapping indices, etc.)
  //
  typedef viennafem::unknown_config<MatrixType,
                                    VectorType,
                                    viennafem::boundary_key,
                                    viennafem::mapping_key>          MyPDEConfigType_1;

  typedef viennafem::unknown_config<MatrixType,
                                    VectorType,
                                    long, char>   MyPDEConfigType_2;  //using keys of type 'long' for boundary, keys of type 'char' for mapping

  MatrixType matrix1, matrix2;
  VectorType rhs1, rhs2;

  MyPDEConfigType_1  poisson_config_1(matrix1,rhs1);
  MyPDEConfigType_2  poisson_config_2(matrix2,rhs2);
  
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
    if (vit->getPoint().get_x() == 0.0 || vit->getPoint().get_x() == 1.0 
      || vit->getPoint().get_y() == 0.0 || vit->getPoint().get_y() == 1.0 )
      viennadata::access<MyPDEConfigType_1::boundary_key_type,
                         bool>(poisson_config_1.boundary_key())(*vit) = true;
    else
      viennadata::access<MyPDEConfigType_1::boundary_key_type,
                         bool>(poisson_config_1.boundary_key())(*vit) = false;
    
    //boundary for second equation: Homogeneous Dirichlet at (x == 0) and (x == 1)
    if (vit->getPoint().get_x() == 0.0 || vit->getPoint().get_x() == 1.0 )
      viennadata::access<MyPDEConfigType_2::boundary_key_type,
                         bool>(poisson_config_2.boundary_key())(*vit) = true;
    else
      viennadata::access<MyPDEConfigType_2::boundary_key_type,
                         bool>(poisson_config_2.boundary_key())(*vit) = false;
  }
  
  
  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_solver fem_solver;

  
  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  fem_solver(poisson_equ_1, poisson_config_1, my_domain);
  fem_solver(poisson_equ_2, poisson_config_2, my_domain);
  
  //std::cout << poisson_config_1.load_vector() << std::endl;
  
  VectorType pde_result_1 = solve(poisson_config_1.system_matrix(), poisson_config_1.load_vector());
  VectorType pde_result_2 = solve(poisson_config_2.system_matrix(), poisson_config_2.load_vector());

  //std::cout << "RESULT" << std::endl;
  //std::cout << pde_result_1 << std::endl;
  //
  // Writing solution back to domain (discussion about proper way of returning a solution required...)
  //
  write_solution_to_VTK_file(pde_result_1, "poisson_1.vtu", my_domain, poisson_config_1);
  write_solution_to_VTK_file(pde_result_2, "poisson_2.vtu", my_domain, poisson_config_2);
  
  std::cout << "*****************************************" << std::endl;
  std::cout << "* Poisson solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
