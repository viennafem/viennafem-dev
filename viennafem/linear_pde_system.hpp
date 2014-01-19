#ifndef VIENNAFEM_LINEAR_PDE_SYSTEM_HPP
#define VIENNAFEM_LINEAR_PDE_SYSTEM_HPP

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

#include <vector>
#include <sstream>

//ViennaFEM includes:
#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"
#include "viennafem/transform.hpp"
#include "viennafem/bases/tetrahedron.hpp"
#include "viennafem/bases/triangle.hpp"
//#include "viennafem/transform/dtdx_triangle.hpp"
//#include "viennafem/transform/dtdx_tetrahedron.hpp"
#include "viennafem/weak_form.hpp"
#include "viennafem/linear_pde_options.hpp"
#include "viennafem/log/latex.hpp"

//ViennaMath includes:
#include "viennamath/manipulation/apply_coordinate_system.hpp"

//ViennaData includes:
#include "viennadata/api.hpp"

//ViennaGrid includes:
#include "viennagrid/config/mesh_config.hpp"


/** @file   linear_pde_system.hpp
    @brief  A class collecting properties of a linear PDE system.
*/

namespace viennafem
{
  /** @brief Representation of a linear system of partial differential equations.
   *
   * @tparam InterfaceType    The ViennaMath runtime expression interface class
   * @tparam MappingKeyType   Type of the ViennaData key used for storing mapping indices
   * @tparam BoundaryKeyType  Type of the ViennaData key used for storing boundary values
   */
  template <typename InterfaceType = viennamath::rt_expression_interface<viennamath::default_numeric_type>,
            typename MappingKeyType  = viennafem::mapping_key,
            typename BoundaryKeyType = viennafem::boundary_key >
  class linear_pde_system
  {
      typedef linear_pde_system<InterfaceType,
                                MappingKeyType,
                                BoundaryKeyType>         self_type;

    public:
      typedef InterfaceType                                 interface_type;
      typedef MappingKeyType                                mapping_key_type;
      typedef BoundaryKeyType                               boundary_key_type;
      typedef viennamath::rt_equation<InterfaceType>           equation_type;
      typedef viennamath::rt_function_symbol<InterfaceType>    unknown_type;
      typedef viennafem::linear_pde_options                 option_type;
      typedef viennafem::latex_logger<InterfaceType>        logger_type;

      linear_pde_system(std::string const & filename) : logger_(filename) {}

      void add_pde(viennamath::rt_equation<InterfaceType> const & pde,
                   std::vector< viennamath::rt_function_symbol<InterfaceType> > const & unknowns,
                   viennafem::linear_pde_options const & option)
      {
        pdes_.push_back(pde);
        unknowns_.push_back(unknowns);
        options_.push_back(option);
      }

      viennamath::rt_equation<InterfaceType> const & pde(std::size_t index) const { return pdes_[index]; }
      std::vector<viennamath::rt_equation<InterfaceType> > const & pdes() const { return pdes_; }

      std::vector<viennamath::rt_function_symbol<InterfaceType> > const & unknown(std::size_t index) const { return unknowns_[index]; }
      viennafem::linear_pde_options option(std::size_t index) const { return options_[index]; }

      std::size_t size() const { return pdes_.size(); }

      logger_type & logger() const { return logger_; }  //TODO: Logger should move to some other place

    private:
      std::vector<viennamath::rt_equation<InterfaceType> >                       pdes_;
      std::vector<std::vector<viennamath::rt_function_symbol<InterfaceType> > >  unknowns_;
      std::vector<viennafem::linear_pde_options>                                 options_;
      mutable logger_type                                                        logger_;
  };

  /** @brief Convenience function for the generation of a linear PDE for a scalar-valued unknown. */
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          viennamath::rt_function_symbol<InterfaceType> unknown_1,
                                                          viennafem::linear_pde_options options_1)
  {
    std::stringstream ss;
    ss << "protocol_" << options_1.data_id() << ".tex";
    linear_pde_system<InterfaceType> ret(ss.str());
    std::vector<viennamath::rt_function_symbol<InterfaceType> > unknown_vec_1(1);
    unknown_vec_1[0] = unknown_1;
    ret.add_pde(equ_1, unknown_vec_1, options_1);
    return ret;
  }

  /** @brief Convenience function for the generation of a linear PDE for a vector-valued unknown. */
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          std::vector<viennamath::rt_function_symbol<InterfaceType> > unknowns_1,
                                                          viennafem::linear_pde_options options_1)
  {
    std::stringstream ss;
    ss << "protocol_" << options_1.data_id() << ".tex";
    linear_pde_system<InterfaceType> ret(ss.str());
    ret.add_pde(equ_1, unknowns_1, options_1);
    return ret;
  }



  /** @brief Convenience function for the generation of a linear PDE for a scalar-valued unknown. */
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          viennamath::rt_function_symbol<InterfaceType> unknown_1)
  {
    linear_pde_system<InterfaceType> ret("protocol_0.tex");

    linear_pde_options options;
    options.data_id(0);
    options.trial_space_id(viennafem::space_to_id<lagrange_tag<1> >::value);
    options.test_space_id(viennafem::space_to_id<lagrange_tag<1> >::value);

    std::vector<viennamath::rt_function_symbol<InterfaceType> > unknown_vec_1(1);
    unknown_vec_1[0] = unknown_1;
    ret.add_pde(equ_1, unknown_vec_1, options);
    return ret;
  }

  /** @brief Convenience function for the generation of a linear PDE for a vector-valued unknown. */
  template <typename InterfaceType>
  linear_pde_system<InterfaceType> make_linear_pde_system(viennamath::rt_equation<InterfaceType> equ_1,
                                                          std::vector<viennamath::rt_function_symbol<InterfaceType> > unknowns_1)
  {
    linear_pde_system<InterfaceType> ret("protocol_0.tex");

    linear_pde_options options;
    options.data_id(0);
    options.trial_space_id(viennafem::space_to_id<lagrange_tag<1> >::value);
    options.test_space_id(viennafem::space_to_id<lagrange_tag<1> >::value);

    ret.add_pde(equ_1, unknowns_1, options);
    return ret;
  }


}
#endif
