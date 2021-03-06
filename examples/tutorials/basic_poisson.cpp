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


#include <iostream>
#include <vector>
#include <stdlib.h>
#include <assert.h>

#include "viennamath/expression.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"

#include "viennamath/runtime/equation.hpp"
#include "viennamath/manipulation/apply_coordinate_system.hpp"

#include "viennafem/weak_form.hpp"

using namespace viennamath;

int main()
{
  std::cout << "************************************************" << std::endl;
  std::cout << "*****     Demo for ViennaMath with FEM     *****" << std::endl;
  std::cout << "************************************************" << std::endl;

  function_symbol u(0, unknown_tag<>());  //an unknown function
  function_symbol v(0, test_tag<>());     //a test function (typical for FEM and FVM)

  //
  // Step 1: Define the classical Poisson equation
  //             laplace(u) = 1
  //
  equation strong_form = make_equation( laplace(u), 1);
  std::cout << "Strong (classical) form of equation:" << std::endl;
  std::cout << "  " << strong_form << std::endl;


  //
  // Step 2: For FEM, an integral formulation ('weak form') is required.
  //         This is obtained by:
  //           - multiplying with a test function (v)
  //           - integration over domain (Omega)
  //           - integration by parts (shift derivative from u to v)
  //
  //         Integral form of Poisson equation:
  //         \int_\Omega grad(u) * grad(v) dx  = \int_\Omega v dx
  //
  equation weak_form_general = viennafem::make_weak_form(strong_form);

  std::cout << "General weak form: " << std::endl;
  std::cout << "  " << weak_form_general << std::endl;

  // Step 3: The final form of the gradient still depends on the underlying coordinate system.
  //         We thus transform the general weak form to cartesian coordinate systems.
  //
  //         Example in 2d:  grad(u) * grad(v)  is transformed to u_x * v_x + u_y * v_y
  //
  equation weak_form_1d = apply_coordinate_system(cartesian<1>(), weak_form_general);
  std::cout << "Weak form in 1d:" << std::endl;
  std::cout << "  " << weak_form_1d << std::endl;

  equation weak_form_2d = apply_coordinate_system(cartesian<2>(), weak_form_general);
  std::cout << "Weak form in 2d:" << std::endl;
  std::cout << "  " << weak_form_2d << std::endl;

  equation weak_form_3d = apply_coordinate_system(cartesian<3>(), weak_form_general);
  std::cout << "Weak form in 3d:" << std::endl;
  std::cout << "  " << weak_form_3d << std::endl;

  //
  // That's it! Everything else is left to ViennaFEM
  //
  std::cout << "************************************************" << std::endl;
  std::cout << "*****     TEST COMPLETED SUCCESSFULLY!     *****" << std::endl;
  std::cout << "************************************************" << std::endl;

  return EXIT_SUCCESS;
}
