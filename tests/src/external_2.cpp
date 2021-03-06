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

//
// A check for the absence of external linkage (otherwise, library is not truly 'header-only')
//

#include <iostream>

//
// *** ViennaFEM
//
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/linear_pde_system.hpp"
#include "viennafem/mapping.hpp"
#include "viennafem/pde_assembler.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/unknown_config.hpp"
#include "viennafem/weak_form.hpp"

#include "viennafem/bases/all.hpp"

#include "viennafem/io/vtk_writer.hpp"

#include "viennafem/log/api.hpp"
#include "viennafem/log/latex.hpp"

#include "viennafem/quadrature/quad.hpp"

#include "viennafem/transform/dtdx_hexahedron.hpp"
#include "viennafem/transform/dtdx_interval.hpp"
#include "viennafem/transform/dtdx_quadrilateral.hpp"
#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/transform/dtdx_triangle.hpp"

void other_func()
{
  viennamath::function_symbol u(0);  //an unknown function

  viennamath::equation strong_form = viennamath::make_equation( viennamath::laplace(u), 1);
  std::cout << "Strong (classical) form of equation:" << std::endl;
  std::cout << "  " << strong_form << std::endl;


  viennamath::equation weak_form_general = viennafem::make_weak_form(strong_form);
}
