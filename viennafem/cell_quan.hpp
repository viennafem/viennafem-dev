/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_CELL_QUAN_HPP
#define VIENNAFEM_CELL_QUAN_HPP

#include "viennafem/forwards.h"

#include "viennamath/forwards.h"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/expression.hpp"
#include "viennadata/api.hpp"

namespace viennafem
{

  template <typename CellType, typename KeyType>
  class cell_quan : public viennamath::default_interface_type
  {
      typedef cell_quan<CellType, KeyType>     self_type;
    public:
      explicit cell_quan(CellType const & cell) : current_cell(&cell) {};
      explicit cell_quan() : current_cell(NULL) {}

      //interface requirements:
      viennamath::default_interface_type * clone() const { return new self_type(*current_cell); }
      viennamath::default_numeric_type eval(std::vector<double> const & v) const
      {
        return viennadata::access<KeyType, viennamath::default_numeric_type>()(*current_cell);
      }
      viennamath::default_numeric_type eval(viennamath::default_numeric_type v) const 
      {
        return viennadata::access<KeyType, viennamath::default_numeric_type>()(*current_cell);
      }
      
      std::string str() const
      {
        std::stringstream ss;
        ss << "cell_quan<" << KeyType() << ">(" << current_cell << ")";
        return ss.str();      
      }
      viennamath::default_numeric_type unwrap() const { throw "Cannot evaluate unknown_func!"; }
      
      viennamath::default_interface_type * substitute(const viennamath::default_interface_type * e,
                                                      const viennamath::default_interface_type * repl) const
      {
        if (dynamic_cast<const self_type *>(e) != NULL)
          return repl->clone();
          
        return clone();
      };    
      
      bool equal(const viennamath::default_interface_type * other) const
      {
        return dynamic_cast< const self_type *>(other) != NULL;
      }
      
      viennamath::default_interface_type * diff(const viennamath::default_interface_type * diff_var) const
      {
        throw "Cannot differentiate cell_quan!";
        return NULL;
      }
      
    private:
      const CellType * current_cell;
  };

  template <typename CellType, typename KeyType>
  viennamath::expr<> operator*(viennamath::variable<> const & lhs,
                             cell_quan<CellType, KeyType> const & rhs)
  {
    return viennamath::expr<>(new viennamath::binary_expr<>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type> >(),
                                                            rhs.clone())); 
  }
  
  
  template <typename CellType, typename KeyType>
  viennamath::expr<> operator*(viennamath::expr<> const & lhs,
                               cell_quan<CellType, KeyType> const & rhs)
  {
    return viennamath::expr<>(new viennamath::binary_expr<>(lhs.get()->clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type> >(),
                                                            rhs.clone())); 
  }
  
  
  
  
  template <typename CellType>
  viennamath::equation<> update_cell_quantities(CellType const & cell,
                                                viennamath::equation<> const & weak_form)
  {
    //step 1: update det_dF_dt
    viennamath::expr<> new_lhs =
       viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType,
                                                                        viennafem::det_dF_dt_key>()),
                              viennamath::expr<>(new viennafem::cell_quan<CellType,
                                                                        viennafem::det_dF_dt_key>(cell)),
                              weak_form.lhs());
                              
    //step 2: update dt_dx<j, i>       
    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 0> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 0> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 1> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 1> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 2> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<0, 2> >(cell)),
                                     new_lhs);
    
    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 0> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 0> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 1> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 1> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 2> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<1, 2> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 0> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 0> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 1> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 1> >(cell)),
                                     new_lhs);

    new_lhs = viennamath::substitute(viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 2> >()),
                                     viennamath::expr<>(new viennafem::cell_quan<CellType, viennafem::dt_dx_key<2, 2> >(cell)),
                                     new_lhs);
    
    
    
    return viennamath::equation<>(new_lhs, weak_form.rhs());
  }
  
  
}
#endif
