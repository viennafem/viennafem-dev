/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_INTERFACE_HPP
#define VIENNAFEM_INTERFACE_HPP

#include "viennafem/pde_assembler.hpp"
#include "viennafem/unknown_config.hpp"
#include "viennamath/expression.hpp"


namespace viennafem
{

template<typename DomainT>
struct assemble
{
   assemble(DomainT& domain) : domain(domain) {}

   template<typename ConfigT>
   void operator()(viennamath::equation<> const& equ, 
                   ConfigT                     & config)
   {
      viennafem::pde_solver assembler;
      std::cout << equ << std::endl;
//      assembler(equ, config, domain);
   }
   DomainT & domain;
};

}


#endif
