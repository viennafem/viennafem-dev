/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                             -----------------

   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaMath base directory
======================================================================= */



#ifndef VIENNAFEM_FEM_EXPRESSION_INTERFACE_HPP
#define VIENNAFEM_FEM_EXPRESSION_INTERFACE_HPP

#include <vector>
#include <string>
#include "viennamath/runtime/binary_expression.hpp"
#include "viennafem/cell_quan.hpp"

namespace viennafem
{
  
  struct cell_updater
  {
    template <typename InterfaceType, typename CellType>
    static void apply(InterfaceType const * p_expr, CellType const & cell)
    {
      if (dynamic_cast<const viennamath::binary_expr<InterfaceType> *>(p_expr) != NULL)
      {
        const viennamath::binary_expr<InterfaceType> * ptr = dynamic_cast<const viennamath::binary_expr<InterfaceType> *>(p_expr);
        ptr->lhs()->update_cell(cell);
        ptr->rhs()->update_cell(cell);
      }
      else if (dynamic_cast<const viennamath::unary_expr<InterfaceType> *>(p_expr) != NULL)
      {
        const viennamath::unary_expr<InterfaceType> * ptr = dynamic_cast<const viennamath::unary_expr<InterfaceType> *>(p_expr);
        ptr->lhs()->update_cell(cell);
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<0,0>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<0,0>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<0,0>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<0,1>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<0,1>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<0,1>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<0,2>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<0,2>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<0,2>" << std::endl;
      }
      

      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<1,0>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<1,0>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<1,0>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<1,1>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<1,1>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<1,1>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<1,2>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<1,2>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<1,2>" << std::endl;
      }
      
      
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<2,0>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<2,0>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<2,0>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<2,1>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<2,1>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<2,1>" << std::endl;
      }
      else if (dynamic_cast<const viennafem::cell_quan<CellType, dt_dx_key<2,2>, InterfaceType> *>(p_expr) != NULL)
      {
        typedef viennafem::cell_quan<CellType, dt_dx_key<2,2>, InterfaceType>    cell_quan_type;
        const cell_quan_type * ptr = dynamic_cast<const cell_quan_type *>(p_expr);
        ptr->update(cell);
        //std::cout << "Updated dt_dx_key<2,2>" << std::endl;
      }
      
      //else
      //  std::cout << "unknown_type" << std::endl;
    }

  };
  
  //class expression_interface;
  //interface for runtime dispatches:
  template <typename NumericT, typename CellType>
  class fem_expression_interface
  {
    public:
      typedef NumericT                                        numeric_type;
      typedef fem_expression_interface<NumericT, CellType>    interface_type;
      
      virtual ~fem_expression_interface() {}
      
      virtual interface_type * clone() const = 0;  //receiver owns pointer!
      virtual std::string str() const = 0;
      virtual NumericT eval(std::vector<NumericT> const & v) const = 0;
      virtual NumericT eval(NumericT val) const = 0;
      virtual interface_type * optimize() const { return clone(); }  //receiver owns pointer!
      virtual bool optimizable() const { return false; }      
      
      virtual bool is_unary() const { return true; }
      
      /** @brief Returns true, if the expression can be evaluated without providing values for variables (i.e. the expression is a constant) */
      virtual bool is_constant() const { return false; }
      
      virtual NumericT unwrap() const = 0;
      virtual interface_type * substitute(const interface_type * e,
                                          const interface_type * repl) const = 0;  //receiver owns pointer! Function parameters must not be manipulated!
                                                
      virtual bool equal(const interface_type * other) const = 0;
      
      virtual const interface_type               * lhs() const { return this; };
      //virtual const op_interface<expression_interface> *          op() const { return NULL; }  //primitives do not have an operator   //TODO: This is a showstopper! Provide something better here!
      //virtual const expression_interface               * rhs() const { return NULL; } //unary expressions do not have a right-hand side


      virtual interface_type * diff(const interface_type * diff_var) const = 0;   //receiver owns pointer! Function parameter diff_var not be manipulated!
      
      
      virtual void update_cell(CellType const & cell) const { cell_updater::apply(this, cell); }
  };
 
}

#endif



