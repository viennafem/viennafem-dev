/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_CELL_QUAN_GUARD
#define VIENNAFEM_CELL_QUAN_GUARD

#include "viennamath/forwards.h"
#include "viennamath/expression.hpp"
#include "viennadata/interface.hpp"

namespace viennafem
{

  template <typename CellType, typename KeyType>
  class cell_quan : public viennamath::expression_interface
  {
      typedef cell_quan<CellType, KeyType>     self_type;
    public:
      explicit cell_quan(CellType const & cell) : current_cell(&cell) {};
      explicit cell_quan() : current_cell(NULL) {}

      //interface requirements:
      expression_interface * clone() const { return new self_type(*current_cell); }
      viennamath::numeric_type eval(std::vector<double> const & v) const
      {
        return viennadata::access<KeyType, viennamath::numeric_type>()(*current_cell);
      }
      viennamath::numeric_type eval(viennamath::numeric_type v) const 
      {
        return viennadata::access<KeyType, viennamath::numeric_type>()(*current_cell);
      }
      
      std::string str() const
      {
        std::stringstream ss;
        ss << "cell_quan<" << KeyType() << ">(" << current_cell << ")";
        return ss.str();      
      }
      viennamath::numeric_type unwrap() const { throw "Cannot evaluate unknown_func!"; }
      
      expression_interface * substitute(const viennamath::expr & e,
                                        const viennamath::expr & repl) const
      {
        if (dynamic_cast<const self_type *>(e.get()) != NULL)
          return repl.get()->clone();
          
        return clone();
      };    
      
      bool equal(const expression_interface * other) const
      {
        return dynamic_cast< const self_type *>(other) != NULL;
      }
      
      expression_interface * diff(const viennamath::expr & diff_var) const
      {
        throw "Cannot differentiate cell_quan!";
        return NULL;
      }
      
    private:
      const CellType * current_cell;
  };

  template <unsigned long id,
            typename CellType, typename KeyType>
  viennamath::expr operator*(viennamath::variable<id> const & lhs,
                             cell_quan<CellType, KeyType> const & rhs)
  {
    return viennamath::expr(new viennamath::binary_expr(lhs.clone(),
                                                        new viennamath::op_mult(),
                                                        rhs.clone())); 
  }
  
  
  template <typename CellType, typename KeyType>
  viennamath::expr operator*(viennamath::expr const & lhs,
                             cell_quan<CellType, KeyType> const & rhs)
  {
    return viennamath::expr(new viennamath::binary_expr(lhs.get()->clone(),
                                                        new viennamath::op_mult(),
                                                        rhs.clone())); 
  }
  
}
#endif
