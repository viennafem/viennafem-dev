#ifndef VIENNAFEM_LOG_INTERFACE_HPP
#define VIENNAFEM_LOG_INTERFACE_HPP

/* =========================================================================
   Copyright (c) 2012, Institute for Microelectronics,
                       Institute for Analysis and Scientific Computing,
                       TU Wien.
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   Author:     Karl Rupp                          rupp@iue.tuwien.ac.at

   License:    MIT (X11), see file LICENSE in the ViennaFEM base directory
============================================================================ */


#include <vector>
#include "viennafem/forwards.h"
#include "viennagrid/forwards.h"
#include "viennagrid/domain.hpp"
#include "viennagrid/segment.hpp"
#include "viennagrid/iterators.hpp"

/** @file  viennafem/log/interface.hpp
    @brief Defines the runtime interface for a logger
*/

namespace viennafem
{
  /** @brief The common interface for all loggers 
   *
   * @tparam InterfaceType   The ViennaMath runtime expression interface
   */
  template <typename InterfaceType>
  class logger_interface
  {
    protected:
      typedef viennamath::rt_equation<InterfaceType>  EquationType;
    
    public:
      virtual void write_strong_form(std::vector<EquationType> const & pdes) = 0;
      virtual void write_weak_form(std::vector<EquationType> const & pdes) = 0;
      virtual void write_coordinated_weak_form(std::vector<EquationType> const & pdes) = 0;
      virtual void write_transformed_weak_form(std::vector<EquationType> const & pdes) = 0;
      virtual void write_test_and_trial_space(std::vector<viennamath::expr> const & test_space,
                                      std::vector<viennamath::expr> const & trial_space) = 0;
      //void write_linear_solver_stats(solver_stats const & solver) = 0;
  };

}

#endif
