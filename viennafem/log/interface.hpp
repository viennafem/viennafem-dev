#ifndef VIENNAFEM_LOG_INTERFACE_HPP
#define VIENNAFEM_LOG_INTERFACE_HPP

/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */


#include <vector>
#include "viennafem/forwards.h"
#include "viennafem/bases/BFStock.hpp"
#include "viennagrid/forwards.h"
#include "viennagrid/domain.hpp"
#include "viennagrid/segment.hpp"
#include "viennagrid/iterators.hpp"

namespace viennafem
{
  
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
