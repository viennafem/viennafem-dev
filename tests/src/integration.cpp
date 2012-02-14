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
#include <iomanip>

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
#include "viennafem/quadrature/quad.hpp"

// ViennaGrid includes:
#include "viennagrid/domain.hpp"
#include "viennagrid/config/simplex.hpp"
#include "viennagrid/config/others.hpp"
#include "viennagrid/io/netgen_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaData includes:
#include "viennadata/api.hpp"

// ViennaMath includes:
#include "viennamath/expression.hpp"



/*
template <typename DomainType>
void fill_with_reference_cell(DomainType & domain, viennagrid::config::line_1d)
{
  typedef viennagrid::config::line_1d  ConfigType;
  
  typedef viennagrid::result_of::point<ConfigType>::type                  PointType;
  typedef viennagrid::result_of::ncell<ConfigType, 0>::type               VertexType;
  typedef viennagrid::result_of::ncell<ConfigType,
                                       ConfigType::cell_tag::dim>::type   CellType;
  
  //push points to domain
  domain.push_back(PointType(0.0));
  domain.push_back(PointType(1.0));
  
  VertexType * cell_vertices[2]; //holds pointers to the respective vertices in the domain
  cell_vertices[0] = &(viennagrid::ncells<0>(domain)[0]);
  cell_vertices[1] = &(viennagrid::ncells<0>(domain)[1]);
  
  CellType cell;
  cell.vertices(cell_vertices);
  domain.push_back(cell);
}

template <typename DomainType>
void fill_with_reference_cell(DomainType & domain, viennagrid::config::triangular_2d)
{
  typedef viennagrid::config::triangular_2d  ConfigType;
  
  typedef viennagrid::result_of::point<ConfigType>::type                  PointType;
  typedef viennagrid::result_of::ncell<ConfigType, 0>::type               VertexType;
  typedef viennagrid::result_of::ncell<ConfigType,
                                       ConfigType::cell_tag::dim>::type   CellType;
  
  //push points to domain
  domain.push_back(PointType(0.0,0.0));
  domain.push_back(PointType(1.0,0.0));
  domain.push_back(PointType(0.0,1.0));
  
  VertexType * cell_vertices[3]; //holds pointers to the respective vertices in the domain
  cell_vertices[0] = &(viennagrid::ncells<0>(domain)[0]);
  cell_vertices[1] = &(viennagrid::ncells<0>(domain)[1]);
  cell_vertices[2] = &(viennagrid::ncells<0>(domain)[2]);
  
  CellType cell;
  cell.vertices(cell_vertices);
  domain.push_back(cell);
}


template <typename DomainType>
void fill_with_reference_cell(DomainType & domain, viennagrid::config::quadrilateral_2d)
{
  typedef viennagrid::config::quadrilateral_2d  ConfigType;
  
  typedef viennagrid::result_of::point<ConfigType>::type                  PointType;
  typedef viennagrid::result_of::ncell<ConfigType, 0>::type               VertexType;
  typedef viennagrid::result_of::ncell<ConfigType,
                                       ConfigType::cell_tag::dim>::type   CellType;
  
  //push points to domain
  domain.push_back(PointType(0.0,0.0));
  domain.push_back(PointType(1.0,0.0));
  domain.push_back(PointType(0.0,1.0));
  domain.push_back(PointType(1.0,1.0));
  
  VertexType * cell_vertices[4]; //holds pointers to the respective vertices in the domain
  cell_vertices[0] = &(viennagrid::ncells<0>(domain)[0]);
  cell_vertices[1] = &(viennagrid::ncells<0>(domain)[1]);
  cell_vertices[2] = &(viennagrid::ncells<0>(domain)[2]);
  cell_vertices[3] = &(viennagrid::ncells<0>(domain)[3]);
  
  CellType cell;
  cell.vertices(cell_vertices);
  domain.push_back(cell);
}

template <typename DomainType>
void fill_with_reference_cell(DomainType & domain, viennagrid::config::hexahedral_3d)
{
  typedef viennagrid::config::hexahedral_3d  ConfigType;
  
  typedef viennagrid::result_of::point<ConfigType>::type                  PointType;
  typedef viennagrid::result_of::ncell<ConfigType, 0>::type               VertexType;
  typedef viennagrid::result_of::ncell<ConfigType,
                                       ConfigType::cell_tag::dim>::type   CellType;
  
  //push points to domain
  domain.push_back(PointType(0.0,0.0,0.0));
  domain.push_back(PointType(1.0,0.0,0.0));
  domain.push_back(PointType(0.0,1.0,0.0));
  domain.push_back(PointType(1.0,1.0,0.0));

  domain.push_back(PointType(0.0,0.0,1.0));
  domain.push_back(PointType(1.0,0.0,1.0));
  domain.push_back(PointType(0.0,1.0,1.0));
  domain.push_back(PointType(1.0,1.0,1.0));
  
  VertexType * cell_vertices[8]; //holds pointers to the respective vertices in the domain
  cell_vertices[0] = &(viennagrid::ncells<0>(domain)[0]);
  cell_vertices[1] = &(viennagrid::ncells<0>(domain)[1]);
  cell_vertices[2] = &(viennagrid::ncells<0>(domain)[2]);
  cell_vertices[3] = &(viennagrid::ncells<0>(domain)[3]);
  
  cell_vertices[4] = &(viennagrid::ncells<0>(domain)[4]);
  cell_vertices[5] = &(viennagrid::ncells<0>(domain)[5]);
  cell_vertices[6] = &(viennagrid::ncells<0>(domain)[6]);
  cell_vertices[7] = &(viennagrid::ncells<0>(domain)[7]);
  
  CellType cell;
  cell.vertices(cell_vertices);
  domain.push_back(cell);
}

template <typename DomainType>
void fill_with_reference_cell(DomainType & domain, viennagrid::config::tetrahedral_3d)
{
  typedef viennagrid::config::tetrahedral_3d  ConfigType;
  
  typedef viennagrid::result_of::point<ConfigType>::type                  PointType;
  typedef viennagrid::result_of::ncell<ConfigType, 0>::type               VertexType;
  typedef viennagrid::result_of::ncell<ConfigType,
                                       ConfigType::cell_tag::dim>::type   CellType;
  
  //push points to domain
  domain.push_back(PointType(0.0,0.0,0.0));
  domain.push_back(PointType(1.0,0.0,0.0));
  domain.push_back(PointType(0.0,1.0,0.0));
  domain.push_back(PointType(0.0,0.0,1.0));
  
  VertexType * cell_vertices[4]; //holds pointers to the respective vertices in the domain
  cell_vertices[0] = &(viennagrid::ncells<0>(domain)[0]);
  cell_vertices[1] = &(viennagrid::ncells<0>(domain)[1]);
  cell_vertices[2] = &(viennagrid::ncells<0>(domain)[2]);
  cell_vertices[3] = &(viennagrid::ncells<0>(domain)[3]);
  
  CellType cell;
  cell.vertices(cell_vertices);
  domain.push_back(cell);
}




template <typename DomainType>
void fill_with_reference_cell(DomainType & d)
{
  typedef typename viennagrid::result_of::config<DomainType>::type   ConfigType;
  
  fill_with_reference_cell(d, ConfigType());
} */


template <typename CellTag>
void fill_integration_rules(std::vector<viennamath::numerical_quadrature> & integrators,
                            std::vector<std::string> & names);

template <>
void fill_integration_rules<viennafem::unit_interval>(std::vector<viennamath::numerical_quadrature> & integrators,
                                                  std::vector<std::string> & names)
{
  typedef viennafem::unit_interval  CellTag;
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 1>()));
  names.push_back("Gauss, exact up to degree 1");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 3>()));
  names.push_back("Gauss, exact up to degree 3");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 5>()));
  names.push_back("Gauss, exact up to degree 5");
  
}


template <>
void fill_integration_rules<viennafem::unit_quadrilateral>(std::vector<viennamath::numerical_quadrature> & integrators,
                                                           std::vector<std::string> & names)
{
  typedef viennafem::unit_quadrilateral  CellTag;
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 1>()));
  names.push_back("Gauss, exact up to degree 1");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 3>()));
  names.push_back("Gauss, exact up to degree 3");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 5>()));
  names.push_back("Gauss, exact up to degree 5");
  
}

template <>
void fill_integration_rules<viennafem::unit_hexahedron>(std::vector<viennamath::numerical_quadrature> & integrators,
                                                        std::vector<std::string> & names)
{
  typedef viennafem::unit_hexahedron  CellTag;
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 1>()));
  names.push_back("Gauss, exact up to degree 1");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 3>()));
  names.push_back("Gauss, exact up to degree 3");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 5>()));
  names.push_back("Gauss, exact up to degree 5");
  
}



template <>
void fill_integration_rules<viennafem::unit_triangle>(std::vector<viennamath::numerical_quadrature> & integrators,
                                                      std::vector<std::string> & names)
{
  typedef viennafem::unit_triangle  CellTag;
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 1>()));
  names.push_back("Gauss, exact up to degree 1");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 7>()));
  names.push_back("Gauss, exact up to degree 7");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 2>()));
  names.push_back("Strang, exact up to degree 2");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 3>()));
  names.push_back("Strang, exact up to degree 3");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 4>()));
  names.push_back("Strang, exact up to degree 4");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 5>()));
  names.push_back("Strang, exact up to degree 5");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 6>()));
  names.push_back("Strang, exact up to degree 6");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 7>()));
  names.push_back("Strang, exact up to degree 7");
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_strang_quad_element<CellTag, 13>()));
  names.push_back("Strang, exact up to degree 13");
  
  
}



template <>
void fill_integration_rules<viennafem::unit_tetrahedron>(std::vector<viennamath::numerical_quadrature> & integrators,
                                                         std::vector<std::string> & names)
{
  typedef viennafem::unit_tetrahedron  CellTag;
  
  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_gauss_quad_element<CellTag, 1>()));
  names.push_back("Gauss, exact up to degree 1");

  //integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 1>()));
  //names.push_back("Keast, exact up to degree 1");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 2>()));
  names.push_back("Keast, exact up to degree 2");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 3>()));
  names.push_back("Keast, exact up to degree 3");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 4>()));
  names.push_back("Keast, exact up to degree 4");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 5>()));
  names.push_back("Keast, exact up to degree 5");

  integrators.push_back(viennamath::numerical_quadrature(new viennafem::rt_keast_quad_element<CellTag, 6>()));
  names.push_back("Keast, exact up to degree 6");

  
}

template <typename RefernceCell>
int test(viennamath::expr const & e)
{
  
  //
  // Run integration:
  //
  std::vector<viennamath::numerical_quadrature> integrators;
  std::vector<std::string> integrator_names;
  
  fill_integration_rules<RefernceCell>(integrators, integrator_names);

  for (std::size_t i=0; i<integrators.size(); ++i)
    std::cout << integrator_names[i] << ": " << std::setprecision(12) << integrators[i](e) << std::endl;

  return EXIT_SUCCESS;
}


int main()
{
  viennamath::variable x(0);
  viennamath::variable y(1);
  viennamath::variable z(2);
  viennamath::constant c(1.0);
  
  /*viennamath::expr e_1d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0) * (x + 2.0) * (x + 3.0) * (x + 4.0));
  viennamath::expr e_2d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0) * (y + 2.0) * (y + 3.0) * (x + 4.0) * (x + 5.0));
  viennamath::expr e_3d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0) * (y + 2.0) * (z + 3.0) + x*y*z);*/
  viennamath::expr e_1d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0));
  viennamath::expr e_2d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0));
  viennamath::expr e_3d = viennamath::integral(viennamath::symbolic_interval(), (x + 1.0));
  
  std::cout << "---- Line ---" << std::endl;
  if (test<viennafem::unit_interval>(e_1d) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::cout << "---- Triangle ---" << std::endl;
  if (test<viennafem::unit_triangle>(e_2d) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::cout << "---- Tetrahedron ---" << std::endl;
  if (test<viennafem::unit_tetrahedron>(e_3d) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::cout << "---- Quadrilateral ---" << std::endl;
  if (test<viennafem::unit_quadrilateral>(e_2d) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::cout << "---- Hexahedron ---" << std::endl;
  if (test<viennafem::unit_hexahedron>(e_3d) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  
  std::cout << "******************************************" << std::endl;
  std::cout << "* Integrator test finished successfully! *" << std::endl;
  std::cout << "******************************************" << std::endl;
 
  return EXIT_SUCCESS;
}