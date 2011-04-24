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

  template <typename CellType, typename KeyType, typename InterfaceType>
  class cell_quan : public InterfaceType
  {
      typedef cell_quan<CellType, KeyType, InterfaceType>     self_type;
      typedef typename InterfaceType::numeric_type            numeric_type;
    public:
      
      explicit cell_quan(CellType const & cell) : current_cell(&cell) {};
      explicit cell_quan() : current_cell(NULL) {}

      //interface requirements:
      InterfaceType * clone() const { return new self_type(*current_cell); }
      numeric_type eval(std::vector<numeric_type> const & v) const
      {
        return viennadata::access<KeyType, numeric_type>()(*current_cell);
      }
      numeric_type eval(numeric_type v) const 
      {
        return viennadata::access<KeyType, numeric_type>()(*current_cell);
      }
      
      std::string str() const
      {
        std::stringstream ss;
        ss << "cell_quan<" << KeyType() << ">(" << current_cell << ")";
        return ss.str();      
      }
      numeric_type unwrap() const { throw "Cannot evaluate unknown_func!"; }
      
      InterfaceType * substitute(const InterfaceType * e,
                                 const InterfaceType * repl) const
      {
        if (dynamic_cast<const self_type *>(e) != NULL)
          return repl->clone();
          
        return clone();
      };    
      
      bool equal(const InterfaceType * other) const
      {
        return dynamic_cast< const self_type *>(other) != NULL;
      }
      
      InterfaceType * diff(const InterfaceType * diff_var) const
      {
        throw "Cannot differentiate cell_quan!";
        return NULL;
      }
      
      
      
      //additional members:
      void update(CellType const & cell) const { current_cell = &cell; }
      
    private:
      mutable const CellType * current_cell;
  };

  template <typename CellType, typename KeyType, typename InterfaceType>
  viennamath::expr<InterfaceType> operator*(viennamath::variable<InterfaceType> const & lhs,
                               cell_quan<CellType, KeyType, InterfaceType> const & rhs)
  {
    return viennamath::expr<InterfaceType>(new viennamath::binary_expr<InterfaceType>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone())); 
  }
  
  
  template <typename CellType, typename KeyType, typename InterfaceType>
  viennamath::expr<InterfaceType> operator*(viennamath::expr<InterfaceType> const & lhs,
                               cell_quan<CellType, KeyType, InterfaceType> const & rhs)
  {
    return viennamath::expr<InterfaceType>(new viennamath::binary_expr<InterfaceType>(lhs.get()->clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone())); 
  }
  
  
  
  /*
  template <typename CellType, typename EquationType>
  viennamath::equation<> update_cell_quantities(CellType const & cell,
                                                EquationType const & weak_form)
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
  } */
  
  
}
#endif
