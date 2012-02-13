/* ====================================================================================
   Copyright (c) 2010, Institute for Microelectronics, Vienna University of Technology.
   http://www.iue.tuwien.ac.at
                                  -----------------
               ViennaFEM - The Vienna Finite Element Method Library
                                  -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaFEM base directory
======================================================================================= */

#ifndef VIENNAFEM_QUADRATURE_TRIANGLE_HPP
#define VIENNAFEM_QUADRATURE_TRIANGLE_HPP

#include "viennafem/forwards.h"
#include "viennafem/cell_quan.hpp"

#include "viennamath/forwards.h"
#include "viennamath/expression.hpp"
#include "viennamath/manipulation/substitute.hpp"
#include "viennamath/manipulation/diff.hpp"
#include "viennamath/manipulation/eval.hpp"
#include "viennamath/runtime/numerical_quadrature.hpp"

#include "viennagrid/topology/triangle.hpp"
#include "viennagrid/topology/tetrahedron.hpp"

namespace viennafem
{
  
  //
  //
  // Exact for polynomials up to order 1
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_triangle, 1, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_triangle, 1, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      explicit rt_gauss_quad_element() : p_(2)
      {
        p_[0] = 1.0/3.0;
        p_[1] = 1.0/3.0;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        return 0.5 * viennamath::eval(e, p_);
      }
      
    private:
      std::vector<NumericT> p_;
  };
  
  
  
  
  //
  //
  // Exact for polynomials up to degree 2
  //
  //
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 2, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 2, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      explicit rt_strang_quad_element() : p1_(2), p2_(2), p3_(2)
      {
        p1_[0] = 2.0/3.0; p1_[1] = 1.0/6.0;
        p2_[0] = 1.0/6.0; p2_[1] = 2.0/3.0;
        p3_[0] = 1.0/6.0; p3_[1] = 1.0/6.0;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        return (viennamath::eval(e, p1_) + viennamath::eval(e, p2_) + viennamath::eval(e, p3_)) / 6.0;
      }
      
    private:
      std::vector<NumericT> p1_;
      std::vector<NumericT> p2_;
      std::vector<NumericT> p3_;
  };
  
  
  
  //
  //
  // Exact for polynomials up to degree 3
  //
  //
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 3, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 3, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 6 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.659027622374092; abscissas_[0][1] = 0.231933368553031;
        abscissas_[1][0] = 0.659027622374092; abscissas_[1][1] = 0.109039009072877;
        abscissas_[2][0] = 0.231933368553031; abscissas_[2][1] = 0.659027622374092;
        abscissas_[3][0] = 0.231933368553031; abscissas_[3][1] = 0.109039009072877;
        abscissas_[4][0] = 0.109039009072877; abscissas_[4][1] = 0.659027622374092;
        abscissas_[5][0] = 0.109039009072877; abscissas_[5][1] = 0.231933368553031;
        
        weights_[0] = 0.16666666666666666667;
        weights_[1] = 0.16666666666666666667;
        weights_[2] = 0.16666666666666666667;
        weights_[3] = 0.16666666666666666667;
        weights_[4] = 0.16666666666666666667;
        weights_[5] = 0.16666666666666666667;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };

  
  
  //
  //
  // Exact for polynomials up to degree 4
  //
  //
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 4, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 4, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 6 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.816847572980459; abscissas_[0][1] = 0.091576213509771;
        abscissas_[1][0] = 0.091576213509771; abscissas_[1][1] = 0.816847572980459;
        abscissas_[2][0] = 0.091576213509771; abscissas_[2][1] = 0.091576213509771;
        abscissas_[3][0] = 0.108103018168070; abscissas_[3][1] = 0.445948490915965;
        abscissas_[4][0] = 0.445948490915965; abscissas_[4][1] = 0.108103018168070;
        abscissas_[5][0] = 0.445948490915965; abscissas_[5][1] = 0.445948490915965;
        
        weights_[0] = 0.109951743655322;
        weights_[1] = 0.109951743655322;
        weights_[2] = 0.109951743655322;
        weights_[3] = 0.223381589678011;
        weights_[4] = 0.223381589678011;
        weights_[5] = 0.223381589678011;
      }
      
      BaseType * clone() const { return new self_type(); }
      
      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  
  
  //
  //
  // Exact for polynomials up to degree 5
  //
  //
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 5, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 5, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 7 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.33333333333333333; abscissas_[0][1] = 0.33333333333333333;
        abscissas_[1][0] = 0.79742698535308720; abscissas_[1][1] = 0.10128650732345633;
        abscissas_[2][0] = 0.10128650732345633; abscissas_[2][1] = 0.79742698535308720;
        abscissas_[3][0] = 0.10128650732345633; abscissas_[3][1] = 0.10128650732345633;
        abscissas_[4][0] = 0.05971587178976981; abscissas_[4][1] = 0.47014206410511505;
        abscissas_[5][0] = 0.47014206410511505; abscissas_[5][1] = 0.05971587178976981;
        abscissas_[6][0] = 0.47014206410511505; abscissas_[6][1] = 0.47014206410511505;
        
        weights_[0] = 0.22500000000000000;
        weights_[1] = 0.12593918054482717;
        weights_[2] = 0.12593918054482717;
        weights_[3] = 0.12593918054482717;
        weights_[4] = 0.13239415278850616;
        weights_[5] = 0.13239415278850616;
        weights_[6] = 0.13239415278850616;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };

  
  //
  //
  // Exact for polynomials up to degree 6
  //
  //
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 6, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 6, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 9 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.124949503233232; abscissas_[0][1] = 0.437525248383384;
        abscissas_[1][0] = 0.437525248383384; abscissas_[1][1] = 0.124949503233232;
        abscissas_[2][0] = 0.437525248383384; abscissas_[2][1] = 0.437525248383384;
        abscissas_[3][0] = 0.797112651860071; abscissas_[3][1] = 0.165409927389841;
        abscissas_[4][0] = 0.797112651860071; abscissas_[4][1] = 0.037477420750088;
        abscissas_[5][0] = 0.165409927389841; abscissas_[5][1] = 0.797112651860071;
        abscissas_[6][0] = 0.165409927389841; abscissas_[6][1] = 0.037477420750088;
        abscissas_[7][0] = 0.037477420750088; abscissas_[7][1] = 0.797112651860071;
        abscissas_[8][0] = 0.037477420750088; abscissas_[8][1] = 0.165409927389841;
        
        weights_[0] = 0.205950504760887;
        weights_[1] = 0.205950504760887;
        weights_[2] = 0.205950504760887;
        weights_[3] = 0.063691414286223;
        weights_[4] = 0.063691414286223;
        weights_[5] = 0.063691414286223;
        weights_[6] = 0.063691414286223;
        weights_[7] = 0.063691414286223;
        weights_[8] = 0.063691414286223;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  
  
  //
  //
  // Exact for polynomials up to degree 7:
  //
  //
  template <typename InterfaceType>
  class rt_gauss_quad_element <viennafem::unit_triangle, 7, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_gauss_quad_element <viennafem::unit_triangle, 7, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 16 };
      
      explicit rt_gauss_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 0.0571041961; abscissas_[0][1] = 0.06546699455602246;
        abscissas_[1][0] = 0.2768430136; abscissas_[1][1] = 0.05021012321401679;
        abscissas_[2][0] = 0.5835904324; abscissas_[2][1] = 0.02891208422223085;
        abscissas_[3][0] = 0.8602401357; abscissas_[3][1] = 0.009703785123906346;
        abscissas_[4][0] = 0.0571041961; abscissas_[4][1] = 0.3111645522491480;
        abscissas_[5][0] = 0.2768430136; abscissas_[5][1] = 0.2386486597440242;
        abscissas_[6][0] = 0.5835904324; abscissas_[6][1] = 0.1374191041243166;
        abscissas_[7][0] = 0.8602401357; abscissas_[7][1] = 0.04612207989200404;
        abscissas_[8][0] = 0.0571041961; abscissas_[8][1] = 0.6317312516508520;
        abscissas_[9][0] = 0.2768430136; abscissas_[9][1] = 0.4845083266559759;
        abscissas_[10][0] = 0.5835904324; abscissas_[10][1] = 0.2789904634756834;
        abscissas_[11][0] = 0.8602401357; abscissas_[11][1] = 0.09363778440799593;
        abscissas_[12][0] = 0.0571041961; abscissas_[12][1] = 0.8774288093439775;
        abscissas_[13][0] = 0.2768430136; abscissas_[13][1] = 0.6729468631859832;
        abscissas_[14][0] = 0.5835904324; abscissas_[14][1] = 0.3874974833777692;
        abscissas_[15][0] = 0.8602401357; abscissas_[15][1] = 0.1300560791760936;
        
        weights_[0] = 0.04713673637581137;
        weights_[1] = 0.07077613579259895;
        weights_[2] = 0.04516809856187617;
        weights_[3] = 0.01084645180365496;
        weights_[4] = 0.08837017702418863;
        weights_[5] = 0.1326884322074010;
        weights_[6] = 0.08467944903812383;
        weights_[7] = 0.02033451909634504;
        weights_[8] = 0.08837017702418863;
        weights_[9] = 0.1326884322074010;
        weights_[10] = 0.08467944903812383;
        weights_[11] = 0.02033451909634504;
        weights_[12] = 0.04713673637581137;
        weights_[13] = 0.07077613579259895;
        weights_[14] = 0.04516809856187617;
        weights_[15] = 0.01084645180365496;
        
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  
  
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 7, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 7, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 13 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[ 0][0] = 0.333333333333333; abscissas_[ 0][1] = 0.333333333333333;
        abscissas_[ 1][0] = 0.479308067841923; abscissas_[ 1][1] = 0.260345966079038;
        abscissas_[ 2][0] = 0.260345966079038; abscissas_[ 2][1] = 0.479308067841923;
        abscissas_[ 3][0] = 0.260345966079038; abscissas_[ 3][1] = 0.260345966079038;
        
        abscissas_[ 4][0] = 0.869739794195568; abscissas_[ 4][1] = 0.065130102902216;
        abscissas_[ 5][0] = 0.065130102902216; abscissas_[ 5][1] = 0.869739794195568;
        abscissas_[ 6][0] = 0.065130102902216; abscissas_[ 6][1] = 0.065130102902216;
        
        abscissas_[ 7][0] = 0.638444188569809; abscissas_[ 7][1] = 0.312865496004875;
        abscissas_[ 8][0] = 0.638444188569809; abscissas_[ 8][1] = 0.048690315425316;
        abscissas_[ 9][0] = 0.312865496004875; abscissas_[ 9][1] = 0.638444188569809;
        abscissas_[10][0] = 0.312865496004875; abscissas_[10][1] = 0.048690315425316;
        abscissas_[11][0] = 0.048690315425316; abscissas_[11][1] = 0.638444188569809;
        abscissas_[12][0] = 0.048690315425316; abscissas_[12][1] = 0.312865496004875;
        
        weights_[ 0] = -0.149570044467670;
        weights_[ 1] = 0.175615257433204;
        weights_[ 2] = 0.175615257433204;
        weights_[ 3] = 0.175615257433204;
        
        weights_[ 4] = 0.053347235608839;
        weights_[ 5] = 0.053347235608839;
        weights_[ 6] = 0.053347235608839;
        
        weights_[ 7] = 0.077113760890257;
        weights_[ 8] = 0.077113760890257;
        weights_[ 9] = 0.077113760890257;
        weights_[10] = 0.077113760890257;
        weights_[11] = 0.077113760890257;
        weights_[12] = 0.077113760890257;
      }

      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  

  
  //
  //
  // Exact for polynomials up to degree 13
  //
  //

  // TOMS algorithm #706
  template <typename InterfaceType>
  class rt_strang_quad_element <viennafem::unit_triangle, 13, InterfaceType> : public viennamath::numerical_quadrature_interface<InterfaceType>
  {
      typedef typename InterfaceType::numeric_type         NumericT;
      typedef rt_strang_quad_element <viennafem::unit_triangle, 13, InterfaceType>  self_type;
      typedef viennamath::numerical_quadrature_interface<InterfaceType>    BaseType;
    public:
      enum { num_points = 37 };
      
      explicit rt_strang_quad_element() : abscissas_(num_points, std::vector<numeric_type>(2)), weights_(num_points)
      {
        abscissas_[0][0] = 1.0/3.0; abscissas_[0][1] = 1.0/3.0;
        
        abscissas_[1][0] = 0.950275662924105565450352089520; abscissas_[1][1] = 0.024862168537947217274823955239;
        abscissas_[2][0] = 0.024862168537947217274823955239; abscissas_[2][1] = 0.950275662924105565450352089520;
        abscissas_[3][0] = 0.024862168537947217274823955239; abscissas_[3][1] = 0.024862168537947217274823955239;
        
        abscissas_[4][0] = 0.171614914923835347556304795551; abscissas_[4][1] = 0.414192542538082326221847602214;
        abscissas_[5][0] = 0.414192542538082326221847602214; abscissas_[5][1] = 0.171614914923835347556304795551;
        abscissas_[6][0] = 0.414192542538082326221847602214; abscissas_[6][1] = 0.414192542538082326221847602214;
        
        abscissas_[7][0] = 0.539412243677190440263092985511; abscissas_[7][1] = 0.230293878161404779868453507244;
        abscissas_[8][0] = 0.230293878161404779868453507244; abscissas_[8][1] = 0.539412243677190440263092985511;
        abscissas_[9][0] = 0.230293878161404779868453507244; abscissas_[9][1] = 0.230293878161404779868453507244;

        abscissas_[10][0] = 0.772160036676532561750285570113; abscissas_[10][1] = 0.113919981661733719124857214943;
        abscissas_[11][0] = 0.113919981661733719124857214943; abscissas_[11][1] = 0.772160036676532561750285570113;
        abscissas_[12][0] = 0.113919981661733719124857214943; abscissas_[12][1] = 0.113919981661733719124857214943;

        abscissas_[13][0] = 0.009085399949835353883572964740; abscissas_[13][1] = 0.495457300025082323058213517632;
        abscissas_[14][0] = 0.495457300025082323058213517632; abscissas_[14][1] = 0.009085399949835353883572964740;
        abscissas_[15][0] = 0.495457300025082323058213517632; abscissas_[15][1] = 0.495457300025082323058213517632;

        abscissas_[16][0] = 0.062277290305886993497083640527; abscissas_[16][1] = 0.468861354847056503251458179727;
        abscissas_[17][0] = 0.468861354847056503251458179727; abscissas_[17][1] = 0.062277290305886993497083640527;
        abscissas_[18][0] = 0.468861354847056503251458179727; abscissas_[18][1] = 0.468861354847056503251458179727;

        abscissas_[19][0] = 0.022076289653624405142446876931; abscissas_[19][1] = 0.851306504174348550389457672223;
        abscissas_[20][0] = 0.022076289653624405142446876931; abscissas_[20][1] = 0.126617206172027096933163647918;
        abscissas_[21][0] = 0.851306504174348550389457672223; abscissas_[21][1] = 0.022076289653624405142446876931;
        abscissas_[22][0] = 0.851306504174348550389457672223; abscissas_[22][1] = 0.126617206172027096933163647918;
        abscissas_[23][0] = 0.126617206172027096933163647918; abscissas_[23][1] = 0.022076289653624405142446876931;
        abscissas_[24][0] = 0.126617206172027096933163647918; abscissas_[24][1] = 0.851306504174348550389457672223;

        abscissas_[25][0] = 0.018620522802520968955913511549; abscissas_[25][1] = 0.689441970728591295496647976487;
        abscissas_[26][0] = 0.018620522802520968955913511549; abscissas_[26][1] = 0.291937506468887771754472382212;
        abscissas_[27][0] = 0.689441970728591295496647976487; abscissas_[27][1] = 0.018620522802520968955913511549;
        abscissas_[28][0] = 0.689441970728591295496647976487; abscissas_[28][1] = 0.291937506468887771754472382212;
        abscissas_[29][0] = 0.291937506468887771754472382212; abscissas_[29][1] = 0.018620522802520968955913511549;
        abscissas_[30][0] = 0.291937506468887771754472382212; abscissas_[30][1] = 0.689441970728591295496647976487;

        abscissas_[31][0] = 0.096506481292159228736516560903; abscissas_[31][1] = 0.635867859433872768286976979827;
        abscissas_[32][0] = 0.096506481292159228736516560903; abscissas_[32][1] = 0.267625659273967961282458816185;
        abscissas_[33][0] = 0.635867859433872768286976979827; abscissas_[33][1] = 0.096506481292159228736516560903;
        abscissas_[34][0] = 0.635867859433872768286976979827; abscissas_[34][1] = 0.267625659273967961282458816185;
        abscissas_[35][0] = 0.267625659273967961282458816185; abscissas_[35][1] = 0.096506481292159228736516560903;
        abscissas_[36][0] = 0.267625659273967961282458816185; abscissas_[36][1] = 0.635867859433872768286976979827;
        
        //weights:
        weights_[0] = 0.051739766065744133555179145422;
        
        weights_[1] = 0.008007799555564801597804123460;
        weights_[2] = 0.008007799555564801597804123460;
        weights_[3] = 0.008007799555564801597804123460;
        
        weights_[4] = 0.046868898981821644823226732071;
        weights_[5] = 0.046868898981821644823226732071;
        weights_[6] = 0.046868898981821644823226732071;
        
        weights_[7] = 0.046590940183976487960361770070;
        weights_[8] = 0.046590940183976487960361770070;
        weights_[9] = 0.046590940183976487960361770070;
        
        weights_[10] = 0.031016943313796381407646220131;
        weights_[11] = 0.031016943313796381407646220131;
        weights_[12] = 0.031016943313796381407646220131;

        weights_[13] = 0.010791612736631273623178240136;
        weights_[14] = 0.010791612736631273623178240136;
        weights_[15] = 0.010791612736631273623178240136;
        
        weights_[16] = 0.032195534242431618819414482205;
        weights_[17] = 0.032195534242431618819414482205;
        weights_[18] = 0.032195534242431618819414482205;
        
        weights_[19] = 0.015445834210701583817692900053;
        weights_[20] = 0.015445834210701583817692900053;
        weights_[21] = 0.015445834210701583817692900053;
        weights_[22] = 0.015445834210701583817692900053;
        weights_[23] = 0.015445834210701583817692900053;
        weights_[24] = 0.015445834210701583817692900053;
        
        weights_[25] = 0.017822989923178661888748319485;
        weights_[26] = 0.017822989923178661888748319485;
        weights_[27] = 0.017822989923178661888748319485;
        weights_[28] = 0.017822989923178661888748319485;
        weights_[29] = 0.017822989923178661888748319485;
        weights_[30] = 0.017822989923178661888748319485;
        
        weights_[31] = 0.037038683681384627918546472190;
        weights_[32] = 0.037038683681384627918546472190;
        weights_[33] = 0.037038683681384627918546472190;
        weights_[34] = 0.037038683681384627918546472190;
        weights_[35] = 0.037038683681384627918546472190;
        weights_[36] = 0.037038683681384627918546472190;
      }
      
      BaseType * clone() const { return new self_type(); }

      NumericT eval(viennamath::rt_interval<InterfaceType> const & interv,
                    viennamath::rt_expr<InterfaceType> const & e,
                    viennamath::rt_variable<InterfaceType> const & var) const
      {
        NumericT result = 0;
        for (std::size_t i=0; i<num_points; ++i)
          result += weights_[i] * viennamath::eval(e, abscissas_[i]);
        return 0.5 * result;
      }
      
    private:
      std::vector<std::vector<NumericT> > abscissas_;
      std::vector<NumericT> weights_;
  };
  
  
}
#endif
