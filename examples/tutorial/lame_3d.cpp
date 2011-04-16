
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>

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
#include "viennagrid/io/sgf_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennamath/runtime/equation.hpp"
#include "viennamath/manipulation/integral.hpp"
#include "viennamath/weak_form.hpp"
#include "viennamath/manipulation/apply_coordinate_system.hpp"


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

using namespace viennamath;

//
// The strain tensor: eps_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
// 
std::vector< expr<> > strain_tensor(std::vector< function_symbol<> > const & u)
{
  //
  // a 3x3 matrix representing the strain tensor
  //
  std::vector< expr<> > result(9);
  
  variable<> x(0);
  variable<> y(1);
  variable<> z(2);
  
  //first row:
  result[0] =        diff(u[0], x);
  result[1] = 0.5 * (diff(u[0], y) + diff(u[1], x));
  result[2] = 0.5 * (diff(u[0], z) + diff(u[2], x));
  
  //second row:
  result[3] = 0.5 * (diff(u[1], x) + diff(u[0], y));
  result[4] =        diff(u[1], y);
  result[5] = 0.5 * (diff(u[1], z) + diff(u[2], y));

  //third row:
  result[6] = 0.5 * (diff(u[2], x) + diff(u[0], z));
  result[7] = 0.5 * (diff(u[2], y) + diff(u[1], z));
  result[8] =        diff(u[2], z);

  return result;
}


//
// The stress tensor: sigma = 2 \mu eps + \lambda trace(eps) Id  for St. Venent-Kirchhoff material
// can be replaced with other expressions for plasticity and the like
// 
std::vector< expr<> > stress_tensor(std::vector< function_symbol<> > const & v)
{
  //
  // a 3x3 matrix representing the stress tensor
  //
  std::vector< expr<> > result(9);
  std::vector< expr<> > strain = strain_tensor(v);

  double mu = 1;
  double lambda = 1;
  
  //The entries are in the following written 
  
  //add 2 \mu eps:
  for (size_t i=0; i<9; ++i)
    result[i] = (2*mu) * strain[i];

  //add trace(eps) * Id:
  result[0] = (2*mu) * strain[0] + lambda * (strain[0] + strain[4] + strain[8]);
  result[4] = (2*mu) * strain[4] + lambda * (strain[0] + strain[4] + strain[8]);
  result[8] = (2*mu) * strain[8] + lambda * (strain[0] + strain[4] + strain[8]);
    
  return result;
}


expr<> tensor_reduce(std::vector< expr<> > lhs, std::vector< expr<> > rhs)
{
  expr<> ret = constant<double>(0);
  
  for (size_t i=0; i<rhs.size(); ++i)
    ret = ret + lhs[i] * rhs[i];
  
  return ret;
}

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
  typedef viennagrid::domain<TetrahedronConfig>         DomainType;

  typedef viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  
  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;
  
  std::cout << "*********************************************************" << std::endl;
  std::cout << "*****     Demo for LAME equation with ViennaFEM     *****" << std::endl;
  std::cout << "*********************************************************" << std::endl;

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

  typedef viennafem::unknown_config<MatrixType,
                                    VectorType,
                                    viennafem::boundary_key,
                                    viennafem::mapping_key>          LameConfig;

  MatrixType matrix;
  VectorType load_vector;

  LameConfig  lame_config(matrix, load_vector);
                                    
                                    
  // the unknown function (vector valued, so one for each of the three components..
  std::vector< function_symbol<> > u(3);
  u[0] = function_symbol<>(0, unknown_tag<>());
  u[1] = function_symbol<>(1, unknown_tag<>());
  u[2] = function_symbol<>(2, unknown_tag<>());
  
  std::vector< function_symbol<> > v(3);
  v[0] = function_symbol<>(0, test_tag<>());
  v[1] = function_symbol<>(1, test_tag<>());
  v[2] = function_symbol<>(2, test_tag<>());
  
  

  //
  // Step 1: Define the classical Lame equation
  //             (lambda + mu) div(u) div(v) + mu grad(u):grad(v) = F
  // with force F set to 0.
  //
  // Minimization problem: \int eps : sigma dx = \int F \cdot u dx
  //

  std::vector< expr<> > strain = strain_tensor(u);
  std::vector< expr<> > stress = stress_tensor(v);
  
  equation<> weak_form_lame = make_equation( 
                                 integral(Omega(), tensor_reduce( strain, stress ), symbolic_tag()),
                                 //=                                         
                                 0); //force is set to zero for now...
  
  
  std::cout << "Weak form of Lame equation: " << std::endl;
  std::cout << weak_form_lame << std::endl;
  
  std::cout << std::endl;
  std::cout << "--- Some bugs in Lame equation handling prohibit further processing. Stop for now... ---" << std::endl;
  std::cout << std::endl;
  exit(0);
  
  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_solver fem_solver;

  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  fem_solver(weak_form_lame, lame_config, my_domain);
  
  VectorType displacements = solve(lame_config.system_matrix(), lame_config.load_vector());

  write_solution_to_VTK_file(displacements, "lame.vtu", my_domain, lame_config);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Lame solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}