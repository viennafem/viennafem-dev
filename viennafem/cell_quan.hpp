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
  
  template <typename CellType, typename DataType>
  class cell_quan_interface
  {
    public: 
      virtual DataType operator()(CellType const & cell) const = 0;
      
      virtual cell_quan_interface<CellType, DataType> * clone() const = 0;
  };

  
  template <typename CellType, typename KeyType, typename DataType>
  class quan_accessor : public cell_quan_interface<CellType, DataType>
  {
    typedef quan_accessor<CellType, KeyType, DataType>    self_type;
    
    public:
      DataType operator()(CellType const & cell) const
      {
        return viennadata::access<KeyType, DataType>()(cell);
      }
      
      cell_quan_interface<CellType, DataType> * clone() const { return new self_type(); }
  };
  
  

  
  template <typename CellType, typename DataType>
  class cell_quan_wrapper
  {
    public:
      template <typename T>
      cell_quan_wrapper(T const * t) : functor_(t) {}
      
      cell_quan_wrapper() {}
      
      cell_quan_wrapper & operator=(cell_quan_wrapper & other)
      {
        functor_ = other.functor_;
        return *this;
      }
      
      DataType operator()(CellType const & cell) const
      {
        return functor_->operator()(cell); 
      }
      
      cell_quan_interface<CellType, DataType> * clone() const { return functor_->clone(); }

    private:
      std::auto_ptr< const cell_quan_interface<CellType, DataType> > functor_;
  };
  
  

  template <typename CellType, typename InterfaceType>
  class cell_quan : public InterfaceType
  {
      typedef cell_quan<CellType, InterfaceType>     self_type;
      typedef typename InterfaceType::numeric_type            numeric_type;
    public:

      explicit cell_quan(CellType const * cell, cell_quan_wrapper<CellType, numeric_type> const & wrapper) : current_cell(cell), accessor(wrapper.clone()) {}
      
      //template <typename T>
      //explicit cell_quan(T const & t) : current_cell(NULL), accessor( new quan_accessor<CellType, T, numeric_type>() ) {}
      
      explicit cell_quan() : current_cell(NULL) {}

      //interface requirements:
      InterfaceType * clone() const { return new self_type(current_cell, accessor); }
      numeric_type eval(std::vector<numeric_type> const & v) const
      {
        return accessor(*current_cell);
      }
      numeric_type eval(numeric_type v) const 
      {
        return accessor(*current_cell);
      }
      
      std::string deep_str() const
      {
        std::stringstream ss;
        ss << "cell_quan(" << current_cell << ")";
        return ss.str();      
      }
      numeric_type unwrap() const { throw "Cannot evaluate unknown_func!"; }
      
      InterfaceType * substitute(const InterfaceType * e,
                                 const InterfaceType * repl) const
      {
        if (deep_equal(e))
          return repl->clone();
          
        return clone();
      };    
      
      InterfaceType * substitute(std::vector<const InterfaceType *> const &  e,
                                 std::vector<const InterfaceType *> const &  repl) const
      {
        //std::cout << "Comparing variable<" << id << "> with " << e->str() << ", result: ";
        for (size_t i=0; i<e.size(); ++i)
          if (deep_equal(e[i]))
            return repl[i]->clone();
        
        //std::cout << "FALSE" << std::endl;
        return clone();
      };    
      
      bool deep_equal(const InterfaceType * other) const
      {
        //TODO: Include comparison of accessor
        return dynamic_cast< const self_type *>(other) != NULL;
      }
      
      bool shallow_equal(const InterfaceType * other) const
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
      
      template <typename T>
      void wrap(T const & t) 
      {
        cell_quan_wrapper<CellType, numeric_type> temp( new quan_accessor<CellType, T, numeric_type>() );
        accessor = temp;
      }
      
    private:
      mutable const CellType * current_cell;
      cell_quan_wrapper<CellType, numeric_type> accessor;
  };

  
  
  template <typename CellType, typename InterfaceType>
  viennamath::expr<InterfaceType> operator*(viennamath::variable<InterfaceType> const & lhs,
                               cell_quan<CellType, InterfaceType> const & rhs)
  {
    return viennamath::expr<InterfaceType>(new viennamath::binary_expr<InterfaceType>(lhs.clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone())); 
  }
  
  
  template <typename CellType, typename InterfaceType>
  viennamath::expr<InterfaceType> operator*(viennamath::expr<InterfaceType> const & lhs,
                               cell_quan<CellType, InterfaceType> const & rhs)
  {
    return viennamath::expr<InterfaceType>(new viennamath::binary_expr<InterfaceType>(lhs.get()->clone(),
                                                            new viennamath::op_binary<viennamath::op_mult<viennamath::default_numeric_type>, InterfaceType >(),
                                                            rhs.clone())); 
  }
  
  
}
#endif
