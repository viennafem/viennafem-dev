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


// remove assert() statements and the like in order to get reasonable performance
#ifndef NDEBUG
  #define NDEBUG
#endif

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>

#include "viennafem/forwards.h"
#include "viennafem/fem.hpp"
#include "viennafem/io/vtk_writer.hpp"

// ViennaGrid includes:
#include "viennagrid/forwards.hpp"
#include "viennagrid/config/default_configs.hpp"
#include "viennagrid/io/netgen_reader.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennamath/runtime/equation.hpp"
#include "viennamath/manipulation/apply_coordinate_system.hpp"


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

#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"

using namespace viennamath;

//
// The strain tensor: eps_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
//
template <typename InterfaceType>
std::vector< rt_expr<InterfaceType> > strain_tensor(std::vector< rt_function_symbol<InterfaceType> > const & u)
{
  typedef rt_variable<InterfaceType>     Variable;

  //
  // a 3x3 matrix representing the strain tensor
  //
  std::vector< rt_expr<InterfaceType> > result(9);

  Variable x(0);
  Variable y(1);
  Variable z(2);

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
template <typename InterfaceType>
std::vector< rt_expr<InterfaceType> > stress_tensor(std::vector< rt_function_symbol<InterfaceType> > const & v)
{
  //
  // a 3x3 matrix representing the stress tensor
  //
  std::vector< rt_expr<InterfaceType> > result(9);
  std::vector< rt_expr<InterfaceType> > strain = strain_tensor(v);

  double mu = 0.5;
  double lambda = 1;

  //The entries are in the following written

  //add 2 \mu eps:
  for (size_t i=0; i<9; ++i)
    result[i] = (2*mu) * strain[i];
    //result[i] = viennamath::constant<>(0);

  //add trace(eps) * Id:
  result[0] = (2*mu) * strain[0] + lambda * (strain[0] + strain[4] + strain[8]);
  result[4] = (2*mu) * strain[4] + lambda * (strain[0] + strain[4] + strain[8]);
  result[8] = (2*mu) * strain[8] + lambda * (strain[0] + strain[4] + strain[8]);

  /*result[0] = lambda * (strain[0] + strain[4] + strain[8]);
  result[4] = lambda * (strain[0] + strain[4] + strain[8]);
  result[8] = lambda * (strain[0] + strain[4] + strain[8]);*/

  return result;
}

//
// Provides the operation a : b, where a and b are tensors
//
template <typename InterfaceType>
rt_expr<InterfaceType> tensor_reduce(std::vector< rt_expr<InterfaceType> > lhs, std::vector< rt_expr<InterfaceType> > rhs)
{
  rt_expr<InterfaceType> ret = lhs[0] * rhs[0];

  for (size_t i=1; i<rhs.size(); ++i)
    ret = ret + lhs[i] * rhs[i];

  return ret;
}


//
// Writes displacements to domain
//
template <typename DomainT, typename StorageT, typename VectorT>
void apply_displacements(DomainT& domain, StorageT& storage, VectorT const & result)
{
  typedef typename viennagrid::result_of::element<DomainT, viennagrid::vertex_tag>::type           VertexType;
  typedef typename viennagrid::result_of::element_range<DomainT, viennagrid::vertex_tag>::type     VertexContainer;
  typedef typename viennagrid::result_of::iterator<VertexContainer>::type                          VertexIterator;

  typedef viennafem::mapping_key          MappingKeyType;
  typedef viennafem::boundary_key         BoundaryKeyType;

  MappingKeyType map_key(0);
  BoundaryKeyType bnd_key(0);

  std::cout << "* apply_displacements(): Writing computed displacements onto domain" << std::endl;
  VertexContainer vertices = viennagrid::elements<VertexType>(domain);
  for (VertexIterator vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    long cur_index = viennadata::access<MappingKeyType, long>(storage, map_key, *vit);
    if (cur_index > -1)
    {
      viennagrid::point(domain, *vit)[0] += result[cur_index+0];
      viennagrid::point(domain, *vit)[1] += result[cur_index+1];
      viennagrid::point(domain, *vit)[2] += result[cur_index+2];
    }
    else
    {
      if (viennadata::access<BoundaryKeyType, std::vector<double> >(storage, bnd_key, *vit).size() > 0)
      {
        viennagrid::point(domain, *vit)[0] += viennadata::access<BoundaryKeyType, std::vector<double> >(storage, bnd_key, *vit)[0];
        viennagrid::point(domain, *vit)[1] += viennadata::access<BoundaryKeyType, std::vector<double> >(storage, bnd_key, *vit)[1];
        viennagrid::point(domain, *vit)[2] += viennadata::access<BoundaryKeyType, std::vector<double> >(storage, bnd_key, *vit)[2];
      }
    }
  }
}

int main()
{
  typedef viennagrid::tetrahedral_3d_mesh                                                 DomainType;
  typedef viennagrid::result_of::segmentation<DomainType>::type                           SegmentationType;
  typedef viennagrid::result_of::element<DomainType, viennagrid::vertex_tag>::type        VertexType;
  typedef viennagrid::result_of::element_range<DomainType, viennagrid::vertex_tag>::type  VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type                          VertexIterator;

  typedef boost::numeric::ublas::compressed_matrix<viennafem::numeric_type>  MatrixType;
  typedef boost::numeric::ublas::vector<viennafem::numeric_type>             VectorType;

  typedef viennamath::function_symbol   FunctionSymbol;
  typedef viennamath::equation          Equation;
  typedef viennamath::expr              Expression;

  typedef viennafem::boundary_key      BoundaryKey;


  std::cout << "*********************************************************" << std::endl;
  std::cout << "*****     Demo for LAME equation with ViennaFEM     *****" << std::endl;
  std::cout << "*********************************************************" << std::endl;

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
    my_reader(my_domain, segments, "../examples/data/cube3072.mesh");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    exit(EXIT_FAILURE);
  }

  MatrixType system_matrix;
  VectorType load_vector;


  // the unknown function (vector valued, so one for each of the three components..
  std::vector< FunctionSymbol > u(3);
  u[0] = FunctionSymbol(0, unknown_tag<>());
  u[1] = FunctionSymbol(1, unknown_tag<>());
  u[2] = FunctionSymbol(2, unknown_tag<>());

  std::vector< FunctionSymbol > v(3);
  v[0] = FunctionSymbol(0, test_tag<>());
  v[1] = FunctionSymbol(1, test_tag<>());
  v[2] = FunctionSymbol(2, test_tag<>());



  //
  // Step 1: Define the classical Lame equation
  //             (lambda + mu) div(u) div(v) + mu grad(u):grad(v) = F
  // with force F set to 0.
  //
  // Minimization problem: \int eps : sigma dx = \int F \cdot u dx
  //

  std::vector< Expression > strain = strain_tensor(u);
  std::vector< Expression > stress = stress_tensor(v);

  Equation weak_form_lame = make_equation( integral(symbolic_interval(), tensor_reduce( strain, stress )),
                                           //=
                                           integral(symbolic_interval(), viennamath::rt_constant<double>(1.0) * v[2])
                                         );


  std::cout << "Weak form of Lame equation: " << std::endl;
  std::cout << weak_form_lame << std::endl;

  std::vector<double> bnd_data_right(3);
  bnd_data_right[0] = 0.2; //small displacement into x-direction prescribed

  VertexContainer vertices = viennagrid::elements<VertexType>(my_domain);
  for (VertexIterator vit = vertices.begin();
      vit != vertices.end();
      ++vit)
  {
    //boundary for first equation: Homogeneous Dirichlet everywhere
    if (viennagrid::point(my_domain, *vit)[0] == 0.0 || viennagrid::point(my_domain, *vit)[0] == 1.0 )
      viennafem::set_dirichlet_boundary(storage, *vit, 0);

    if (viennagrid::point(my_domain, *vit)[0] == 1.0)
    {
      viennafem::set_dirichlet_boundary(storage, *vit, bnd_data_right);
      viennadata::access<BoundaryKey, double>(storage, BoundaryKey(0), *vit) = bnd_data_right[0]; //this is for the moment used for the VTK writer
    }
  }

  //
  // Create PDE solver functors: (discussion about proper interface required)
  //
  viennafem::pde_assembler<StorageType> fem_assembler(storage);

  //
  // Assemble and solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  fem_assembler(viennafem::make_linear_pde_system(weak_form_lame,
                                                  u,
                                                  viennafem::make_linear_pde_options(0,
                                                                                     viennafem::lagrange_tag<1>(),
                                                                                     viennafem::lagrange_tag<1>())
                                                 ),
                my_domain,
                system_matrix,
                load_vector
               );

  VectorType displacements = viennacl::linalg::solve(system_matrix, load_vector, viennacl::linalg::bicgstab_tag());
  std::cout << "* solve(): Residual: " << norm_2(prod(system_matrix, displacements) - load_vector) << std::endl;

  apply_displacements(my_domain, storage, displacements);
  viennafem::io::write_solution_to_VTK_file(displacements, "lame", my_domain, segments, storage, 0);

  std::cout << "*****************************************" << std::endl;
  std::cout << "* Lame solver finished successfully! *" << std::endl;
  std::cout << "*****************************************" << std::endl;
  return EXIT_SUCCESS;
}
