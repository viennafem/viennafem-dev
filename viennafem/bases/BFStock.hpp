/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_BFSTOCK_HPP
#define VIENNAFEM_BFSTOCK_HPP

#include <vector>
#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"
#include "viennagrid/forwards.h"
#include "viennamath/expression.hpp"

namespace viennafem
{

  
  
  
  
  /////////////////// Old code to follow below this line ////////////////////////////////////////
  
  
  
  
  
  
  

  //helper metafunction (that also offers runtime-version) for determination of degrees of freedom
  template <long degree = 0,
            long tuple_len = 0,
            bool args_illegal = (degree < 0) || (tuple_len < 0) || (degree - tuple_len < 0) >
  class PascalSimplex
  {
    public:
      enum{ ReturnValue = PascalSimplex<degree - 1, tuple_len - 1>::ReturnValue
                          +PascalSimplex<degree - 1, tuple_len>::ReturnValue };
  };

  template <long degree>
  class PascalSimplex<degree, degree, false>
  {
    public:
      enum{ ReturnValue = 1 };
  };

  template <long degree>
  class PascalSimplex<degree, 1, false>
  {
    public:
      enum{ ReturnValue = 1 };
  };

  template <>
  class PascalSimplex<1, 1, false>
  {
    public:
      enum{ ReturnValue = 1 };
  };

  template <>
  class PascalSimplex<0,0, false>
  {
    public:

      enum{ ReturnValue = 0 };

      static PascalSimplex & getInstance()
      {
        static PascalSimplex ps;
        return ps;
      }

      long getNum(long degree_, long tuple_len) const
      {
        if (tuple_len == 1)
          return 1;
        if (tuple_len == degree_)
          return 1;

        return getNum(degree_-1, tuple_len-1) + getNum(degree_-1, tuple_len);
      }

    private:
      PascalSimplex() {};
      PascalSimplex( const PascalSimplex &) {};
  };

  //prevent potential abuse via negative arguments:
  template <long degree, long tuple_len>
  class PascalSimplex<degree, tuple_len, true>
  {
    public:
      enum{ ReturnValue = 0 };
  };

  
  
  template <typename CellTag>
  struct BFStock
  {};

  /*
  template <>
  struct BFStock<viennagrid::point_tag>
  {
      //provide trivial compound type:
      typedef ScalarExpression<1>                          CompoundType;

  };

  template <>
  struct BFStock<LineTag>
  {
      //provide CompoundType:
      typedef Expression < ExpressionDefaultScalarType,
                            ScalarExpression<1>,
                            var<0>,
                            op_minus<ExpressionDefaultScalarType> >                    OneMinusX;
      typedef CompoundExpression<1, 1, OneMinusX, 0, var<0>, 0>     CompoundType;

  };

  template <>
  struct BFStock<TriangleTag>
  {
      //provide CompoundType:
      typedef Expression < ExpressionDefaultScalarType,
                            ScalarExpression<1>,
                            Expression< ExpressionDefaultScalarType,
                                        var<0>,
                                        var<1>,
                                        op_plus<ExpressionDefaultScalarType>
                                      >,
                            op_minus<ExpressionDefaultScalarType> >                    OneMinusXY;
      typedef CompoundExpression<1, 1, OneMinusXY, 0, var<0>, 0, var<1>, 0>     CompoundType;

  }; //BFStock<triangle>


  template <>
  struct BFStock<TetrahedronTag>
  {
      typedef Expression < ExpressionDefaultScalarType,
                            Expression< ExpressionDefaultScalarType,
                                        ScalarExpression<1>,
                                        var<0>,
                                        op_minus<ExpressionDefaultScalarType>
                                      >,
                            Expression< ExpressionDefaultScalarType,
                                        var<1>,
                                        var<2>,
                                        op_plus<ExpressionDefaultScalarType>
                                      >,
                            op_minus<ExpressionDefaultScalarType> >                      OneMinusXYZ;
      typedef CompoundExpression<1, 1, OneMinusXYZ, 0, var<0>, 0, var<1>, 0, var<2>, 0>  CompoundType;

  }; */

  template <long topolevel_, long element_id_, long bf_id_>
  struct BasisFunctionID
  {
    enum{ topolevel = topolevel_,
           element_id = element_id_,
           bf_id = bf_id_ };

    void print()
    {
      std::cout << "<" << topolevel << "," << element_id << "," << bf_id << ">" << std::endl;
    }
  };

  struct BasisFunctionIDInvalid {};   //end of iteration
  struct IteratorStop {};     //stopping criterion for compiletime iterators:

  /*
  template <typename CellTag, long element_id, long bf_id,
              typename ElementTag, long bf_degree>
  struct GET_COMPRESSED_BF
  {
    //build compressed_basisfun for element's vertex-functions:
    typedef typename GET_COMPRESSED_BF_FOR_ELEMENT
                                          < TopologyLevel<ElementTag, 0>::ElementNum,
                                            bf_id,
                                            bf_degree>::ResultType            ElementBfType;


    //embed this element-compressed_basisfun into cell-compressed_basisfun:
    typedef typename EMBED_ELEMENT_BF_TO_CELL_BF<ElementTag,
                                                    CellTag,
                                                    ElementBfType,
                                                    element_id
                                                    >::ResultType   ResultType;

  };

  //all-in-one interface: returns the desired basisfunction:
  template <typename AssemblyCellTag, typename BaseFunTag,
              typename BFID>
  struct GET_BASISFUNCTION
  {
    //construct the desired compressed_basisfun:
    typedef typename GET_COMPRESSED_BF<AssemblyCellTag,
                                          BFID::element_id,
                                          BFID::bf_id,
                                          typename TopologyLevel<AssemblyCellTag,
                                                                  BFID::topolevel>::ElementTag,
                                          BaseFunTag::degree
                                         >::ResultType                  CompressedBfType;

    //finally, the compressed_basisfun is transformed to a valid expression:
    typedef typename COMPRESSED_BF_TO_FULL_BF<CompressedBfType, AssemblyCellTag>::ResultType       ResultType;

  }; */


  ///////////////// Support for type erasure /////////////////

  //forward declaration of helper metafunction:
  template <long element_id, long bf_id, long bf_per_element >
  struct MatrixIteratorIncrementer;

  //Helper class for conversion basisfunction-typelist -> basisfunction-array
  template <typename FEMConfig, typename AssemblyCellType,
              typename BFIterator, long bf_array_index = 0
              >
  struct BFFiller
  {
    template <typename BasisFunction>
    static BasisFunction * apply(BasisFunction * basisfuns)
    {
      typedef typename BFIterator::ResultType   BFType;

      BFType expr;
      BasisFunction bf( expr );
      basisfuns[bf_array_index] = bf;
      BFFiller< FEMConfig, AssemblyCellType,
                typename BFIterator::NextType,
                bf_array_index + 1
              >::apply(basisfuns);
      return basisfuns;
    }
  };

  template <typename FEMConfig, typename AssemblyCellType, long bf_array_index>
  struct BFFiller<AssemblyCellType, FEMConfig, IteratorStop, bf_array_index>
  {
    template <typename BasisFunction>
    static void apply(BasisFunction * basisfuns) {}
  };


  // Basisfunction accessor:
  // holds index-tuple (i1, i2,...,ik), where sum(i_j) = d
  // can be seen as a std::vector with reduced functionality, but additional operator++
  class BaseFunKey
  {
    public:
      BaseFunKey(short tuple_length, short bf_degree) 
        : key(tuple_length), degree(bf_degree), valid_(true)
      {
        //set first key:
        //std::cout << "Constructing BaseFunKey" << std::endl;
        for (short i=0; i<degree; ++i)
          key[i] = 1;
        key[tuple_length-1] = bf_degree - tuple_length + 1;
      };

      BaseFunKey(const BaseFunKey & bfk) : key(bfk.key.size()), degree(bfk.degree), valid_(bfk.valid_)
      {
        for (short i=0; i<degree; ++i)
          key[i] = bfk.key[i];
      }


      long operator[](unsigned long index) const { return key[index]; }

      BaseFunKey & operator++()
      {
        short temp = 0;
        short inc_bound = degree - key.size() + 1;

        //issue warning if not valid anymore:
        if (!valid_)
          std::cerr << "BaseFunKey operator++: end already reached" << std::endl;

        //find next increment-decrement-pair:
        for (short i=degree - 1; i > 0; --i)
        {
          if ( (key[i] > 1) && (key[i-1] < inc_bound) )
          {
            key[i-1] += 1;
            temp = key[i]; //temp avoids checking for (i == degree_ - 1)
            key[i] = 1;

            //balance last index:
            key[key.size() - 1] = temp - 1;

            return *this;
          }
        }

        valid_ = false;
        return *this;
      }

      //no operator++(int), since array-copying is too costly

      void print()
      {
        for (unsigned short i=0; i<key.size(); ++i)
          std::cout << key[i] << ",";
        std::cout << std::endl;
      }

      short getDegree() const { return degree; }

      short size() const { return key.size(); }

      bool valid() const { return valid_; }

    private:
      std::vector<short> key;
      short degree;
      bool valid_;
  };



  ///////// Iterator for basisfunction IDs /////////////////////

  //incrementer for bf_id:
  template <typename AssemblyCellType, typename BasisfunctionTag, 
              long topolevel, long current_bf_id,
              bool increment_possible = true
                  //(BasisFuncNum< typename AssemblyCellType::ElementTag,
                  //               BasisfunctionTag,
                  //               topolevel>::ReturnValue - current_bf_id > 1)
            >
  struct BfIDIncrementer_bf_id
  {
    enum { Incrementable = 1 };
  };

  template <typename AssemblyCellType, typename BasisfunctionTag, 
              long topolevel, long current_bf_id>
  struct BfIDIncrementer_bf_id<AssemblyCellType, BasisfunctionTag, topolevel, current_bf_id, false>
  {
    enum { Incrementable = 0 };
  };

  //incrementer for element_id:
  template <typename AssemblyCellType,
              long topolevel, long current_element_id,
              bool increment_possible =
                  (viennagrid::topology::subcell_desc< typename AssemblyCellType::ElementTag,
                                  topolevel
                                >::num_elements - current_element_id > 1)
            >
  struct BfIDIncrementer_element_id
  {
    enum { Incrementable = 1 };
  };

  template <typename AssemblyCellType,
              long topolevel, long current_element_id>
  struct BfIDIncrementer_element_id<AssemblyCellType, topolevel, current_element_id, false>
  {
    enum { Incrementable = 0 };
  };

  //incrementer for topolevel:
  template <typename AssemblyCellType, typename BasisfunctionTag, 
              long topolevel,
              bool next_topolevel_possible = true,
//                  (BasisFuncNum< typename AssemblyCellType::ElementTag,
//                                 BasisfunctionTag,
//                                topolevel + 1>::ReturnValue > 0),
              bool topolevel_valid = (AssemblyCellType::ElementTag::TopoLevel - topolevel > 0)
            >
  struct BfIDIncrementer_topolevel
  {
    //standard case: increment topolevel by one.
    enum { Incrementable = 1,
            Next_topolevel = topolevel + 1 };
  };

  template <typename AssemblyCellType, typename BasisfunctionTag, 
              long topolevel>
  struct BfIDIncrementer_topolevel<AssemblyCellType, BasisfunctionTag, topolevel, false, true>
  {
    //try next level:
    enum { Incrementable = BfIDIncrementer_topolevel<AssemblyCellType,
                                                      BasisfunctionTag,
                                                      topolevel + 1>::Incrementable,
            Next_topolevel = BfIDIncrementer_topolevel<AssemblyCellType,
                                                      BasisfunctionTag,
                                                      topolevel + 1>::Next_topolevel };
  };

  template <typename AssemblyCellType, typename BasisfunctionTag, 
              long topolevel, bool next_topolevel_possible>
  struct BfIDIncrementer_topolevel<AssemblyCellType, BasisfunctionTag, topolevel,
                                    next_topolevel_possible, false>
  {
    //end recursive topolevel-search:
    enum { Incrementable = 0, Next_topolevel = topolevel };
  };


  template <typename AssemblyCellType, typename BasisfunctionTag, typename BfID,
              bool inc_topolevel_possible =
                    (BfIDIncrementer_topolevel<AssemblyCellType, BasisfunctionTag, 
                                          BfID::topolevel>::Incrementable > 0),
              bool inc_element_possible =
                    (BfIDIncrementer_element_id<AssemblyCellType,
                                          BfID::topolevel, BfID::element_id>::Incrementable > 0),
              bool inc_bf_possible =
                    (BfIDIncrementer_bf_id<AssemblyCellType, BasisfunctionTag, 
                                          BfID::topolevel, BfID::bf_id>::Incrementable > 0)
            >
  struct BfIDIncrementer
  {
    //default case: increment bf:
    typedef BasisFunctionID<BfID::topolevel,
                             BfID::element_id,
                             BfID::bf_id + 1>     ResultType;
  };

  //increment element_id (when bf_id cannot be incremented anymore)
  template <typename AssemblyCellType, typename BasisfunctionTag, typename BfID,
              bool inc_topolevel_possible>
  struct BfIDIncrementer<AssemblyCellType, BasisfunctionTag, BfID,
                          inc_topolevel_possible, true, false>
  {
    typedef BasisFunctionID<BfID::topolevel,
                             BfID::element_id + 1,
                             0>     ResultType;
  };

  //increment topolevel (when element_id and bf_id cannot be incremented anymore)
  template <typename AssemblyCellType, typename BasisfunctionTag, typename BfID>
  struct BfIDIncrementer<AssemblyCellType, BasisfunctionTag, BfID,
                          true, false, false>
  {
    typedef BasisFunctionID<BfIDIncrementer_topolevel<AssemblyCellType, BasisfunctionTag, 
                                                        BfID::topolevel>::Next_topolevel,
                             0,
                             0>     ResultType;
  };

  //stop iteration if no more increments possible:
  template <typename AssemblyCellType, typename BasisfunctionTag, typename BfID>
  struct BfIDIncrementer<AssemblyCellType, BasisfunctionTag, BfID,
                          false, false, false>
  {
    typedef BasisFunctionIDInvalid     ResultType;
  };
  

  template <typename AssemblyCellType, typename BasisfunctionTag, 
              typename BfID = BasisFunctionID<0,0,0> >
  struct BasisFunctionIterator;

  //terminator:
  template <typename AssemblyCellType, typename BasisfunctionTag, typename NextBFID>
  struct BFIteratorNext
  {
    typedef BasisFunctionIterator<AssemblyCellType, BasisfunctionTag, NextBFID>   ResultType;
  };

  template <typename AssemblyCellType, typename BasisfunctionTag>
  struct BFIteratorNext<AssemblyCellType, BasisfunctionTag, BasisFunctionIDInvalid>
  {
    typedef IteratorStop   ResultType;
  };

  //All-In-One basisfunction-ID iterator
  template <typename AssemblyCellType, typename BfTag, typename BfID>
  struct BasisFunctionIterator
  {
    typedef BfTag                   BasisfunctionTag;
    typedef BfID                    CurrentID;

//    typedef typename GET_BASISFUNCTION<typename AssemblyCellType::ElementTag,
//                                          BasisfunctionTag,
//                                         CurrentID>::ResultType   ResultType;

    typedef typename BFIteratorNext<AssemblyCellType, BasisfunctionTag, 
                   typename BfIDIncrementer<AssemblyCellType, BasisfunctionTag, BfID>::ResultType
                                      >::ResultType                         NextType;
  };


}

#endif
