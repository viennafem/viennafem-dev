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
#include "viennafem/pde_config.hpp"
#include "viennafem/pde_solver.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include "viennagrid/io/sgf_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/interface.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"
#include "viennamath/equation.hpp"
#include "viennamath/unknown_func.hpp"
#include "viennamath/op_tags.hpp"


//
// Configuration class for a triangular domain
//
struct TriangleConfig
{
  typedef double                                  numeric_type;
  typedef viennagrid::two_dimensions_tag          dimension_tag;
  typedef viennagrid::triangle_tag                cell_tag;

  //multigrid:
  //typedef viennagrid::full_multigrid_tag                       multigrid_tag;
  typedef viennagrid::no_multigrid_tag             multigrid_tag;
};


template <typename DomainType,
          typename PDEConfig>
void write_solution_to_VTK_file(std::vector<viennafem::numeric_type> const & result,
                                std::string filename,
                                DomainType & domain,
                                PDEConfig const & config)
{
  typedef typename viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
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
  typedef viennagrid::domain<TriangleConfig>         DomainType;

  typedef viennagrid::result_of::ncell_container<DomainType, 0>::type    VertexContainer;
  typedef viennagrid::result_of::iterator<VertexContainer>::type         VertexIterator;
  
  
  //
  // Create a domain from file
  //
  DomainType my_domain;
  
  try
  {
    viennagrid::io::sgf_reader my_sgf_reader;
    my_sgf_reader(my_domain, "../examples/data/square128.sgf");
  }
  catch (...)
  {
    std::cerr << "File-Reader failed. Aborting program..." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  
  //
  // Specify two PDEs:
  //
  viennamath::unknown_func<0> u;   //an unknown function used for PDE specification
  viennamath::equation poisson_equ_1 = viennamath::make_equation( viennamath::laplace(u), -1);
  viennamath::equation poisson_equ_2 = viennamath::make_equation( viennamath::laplace(u), -1);
  
  //
  // Create PDE config (where to find boundary information, where to store mapping indices, etc.)
  //
  typedef viennafem::pde_config<viennafem::boundary_key,
                                viennafem::mapping_key>          MyPDEConfigType_1;

  typedef viennafem::pde_config<long, char>   MyPDEConfigType_2;  //using keys of type 'long' for boundary, keys of type 'char' for mapping
                                
  MyPDEConfigType_1  poisson_config_1;
  MyPDEConfigType_2  poisson_config_2;
  
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
  viennafem::pde_solver<MyPDEConfigType_1> my_solver_1(poisson_config_1);
  viennafem::pde_solver<MyPDEConfigType_2> my_solver_2(poisson_config_2);

  //
  // Solve system and write solution vector to pde_result:
  // (discussion about proper interface required. Introduce a pde_result class?)
  //
  std::vector<viennafem::numeric_type> pde_result_1 = my_solver_1(poisson_equ_1, my_domain);
  std::vector<viennafem::numeric_type> pde_result_2 = my_solver_2(poisson_equ_2, my_domain);
  
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
