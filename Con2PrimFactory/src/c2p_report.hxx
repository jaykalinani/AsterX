// Error report adapted from the RePrimAnd library

#ifndef C2P_REPORT_HXX
#define C2P_REPORT_HXX

namespace Con2PrimFactory {

//class con2prim_mhd;
//namespace detail {class froot;}

///Struct to represent the outcome of a con2prim call.
/**
This should be checked after con2prim. The various members indicate
success or failure, if the conserved variables were changed, or
artificial atmosphere applied.
*/
struct c2p_mhd_report  {
  public:
  
  ///Failure flags 
  enum err_code {
    SUCCESS,       /*!<Recovery succeeded (possibly after 
                       allowed corrections) */
    INVALID_DETG,  ///<3-metric determinant not positive or NAN/INF
    NEG_BSQR,      ///<Negative B square, metric not positive definite 
    NANS_IN_CONS,  ///<One or more evolved variables is NAN/INF
    RANGE_RHO,     ///<Mass density outside EOS range
    RANGE_EPS,     ///<Fluid internal energy outside EOS range
    SPEED_LIMIT,   ///<Speed limit exceeded
    RANGE_YE,      ///<Electron fraction outside EOS range
    B_LIMIT,       ///<Magnetization above limit
    ROOT_FAIL_CONV,        /*!<Root finding did not converge (this 
                               should never happen) */ 
    ROOT_FAIL_BRACKET,     /*!<Bracketing of root failed (this 
                               should never happen) */ 
    PREP_ROOT_FAIL_CONV,   /*!<Auxiliary root finding did not 
                               converge (this should never happen) */ 
    PREP_ROOT_FAIL_BRACKET,/*!<Bracketing of auxiliary root 
                               failed (this should never happen) */ 
    ERR_CODE_NOT_SET
  };
  
  ///Default constructor, resulting object invalid.
  c2p_mhd_report() = default;
  
  
  /**
  @return String with debug information.
  
  This assembles a human readable debug message from the flags
  and numerical values stored internally (for performance reasons, no
  string is ever created before explicitly requested via this function).
  **/
  void debug_message() const;
  
  /** 
  @return If the input was invalid according to the error policy.
  **/
  bool failed() const {return status != SUCCESS;}        
  
  
  /// SUCCESS or reason for failure. 
  err_code status{ERR_CODE_NOT_SET};

  /// Whether the conserved variables were adjusted.
  bool adjust_cons{false};     

  /// Whether the artificial atmosphere was enforced.  
  bool set_atmo{false};        

  /// Number of calls to the EOS needed for the root finding.
  unsigned CCTK_INT iters{0};            

  protected:

  /// Set state artificial atmosphere enforced.
  void set_atmo_set();
  
  /// Set error invalid metric determinant. 
  void set_invalid_detg(CCTK_REAL detg_);
  
  /// Set error negative B square. 
  void set_neg_bsqr(CCTK_REAL bsqr_);
  
  /// Set error NANs in conserved variables. 
  void set_nans_in_cons(CCTK_REAL dens_, CCTK_REAL qtot_, CCTK_REAL rsqr_,
    CCTK_REAL rbsqr_, CCTK_REAL bsqr_, CCTK_REAL ye_);
    
  /// Set error density out of range. 
  void set_range_rho(CCTK_REAL dens_, CCTK_REAL rho_);
    
  /// Set error energy out of range. 
  void set_range_eps(CCTK_REAL eps_);
  
  /// Set error speed limit exceeded. 
  void set_speed_limit(CCTK_REAL vel_);

  /// Set error limit for b exceeded. 
  void set_b_limit(CCTK_REAL bsqr_);
  
  /// Set error energy out of range. 
  void set_range_ye(CCTK_REAL ye_);
  
  /// Set error root finding did not converge. 
  void set_root_conv();

  /// Set error root bracketing faulty. 
  void set_root_bracket();
  
  /// Set error preperatory root finding not converged. 
  void set_prep_root_conv();

  /// Set error preperatory root bracketing faulty. 
  void set_prep_root_bracket();

  private:
  
  /// Conserved density \f$ D \f$.
  CCTK_REAL dens;

  /// Specific conserved momentum \f$ r^2 = \frac{ S_i S^i}{D^2} \f$.
  CCTK_REAL rsqr;               

  /// Product \f$ (r^l b_l)^2 \f$.
  CCTK_REAL rbsqr;              

  /// Specific total energy density \f$ q = \frac{\tau}{D}  \f$.
  CCTK_REAL qtot;               

  /// Ratio \f$ b^2 = \frac{B^2}{D} \f$.
  CCTK_REAL bsqr;

  /// Electron fraction \f$ Y_e \f$.
  CCTK_REAL ye;
  
  /// Rest mass density \f$ \rho \f$.
  CCTK_REAL rho;
  
  /// Specific internal energy \f$ \epsilon \f$.
  CCTK_REAL eps;
  
  /// Velocity.
  CCTK_REAL vel;
  
  ///3-Metric determinant.
  CCTK_REAL detg; 
  
  
//  friend class con2prim_mhd;
//  friend class detail::froot;

// CXX implementation

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_atmo_set()
{ 
  status      = SUCCESS;
  set_atmo    = true;
  adjust_cons = true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_invalid_detg(real_t detg_)
{ 
  status      = INVALID_DETG;
  adjust_cons = true;
  detg        = detg_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_neg_bsqr(real_t bsqr_)
{ 
  status      = NEG_BSQR;
  adjust_cons = true;
  bsqr        = bsqr_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_nans_in_cons(real_t dens_, real_t qtot_,
  real_t rsqr_, real_t rbsqr_, real_t bsqr_, real_t ye_)
{ 
  status      = NANS_IN_CONS;
  set_atmo    = false;
  adjust_cons = true;
  dens        = dens_;
  qtot        = qtot_;
  rsqr        = rsqr_;
  rbsqr       = rbsqr_;
  bsqr        = bsqr_;
  ye          = ye_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_range_rho(real_t dens_, real_t rho_)
{ 
  status      = RANGE_RHO;
  set_atmo    = false;
  adjust_cons = true;
  dens        = dens_;
  rho         = rho_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_range_eps(real_t eps_)
{
  status      = RANGE_EPS;
  set_atmo    = false;
  adjust_cons = true;
  eps         = eps_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_speed_limit(real_t vel_)
{
  status      = SPEED_LIMIT;
  set_atmo    = false;
  adjust_cons = true;
  vel         = vel_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_b_limit(real_t bsqr_)
{
  status      = B_LIMIT;
  set_atmo    = false;
  adjust_cons = true;
  bsqr        = bsqr_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_range_ye(real_t ye_)
{
  status      = RANGE_YE;
  set_atmo    = false;
  adjust_cons = true;
  ye          = ye_;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_root_conv()
{
  status      = ROOT_FAIL_CONV;
  set_atmo    = false;
  adjust_cons = true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_root_bracket()
{
  status      = ROOT_FAIL_BRACKET;
  set_atmo    = false;
  adjust_cons = true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_prep_root_conv()
{
  status          = PREP_ROOT_FAIL_CONV;
  set_atmo        = false;
  adjust_cons     = true;
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_prep_root_bracket()
{
  status          = PREP_ROOT_FAIL_BRACKET;
  set_atmo        = false;
  adjust_cons     = true;
}


CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void debug_message() const
{
  switch (status) {
    case SUCCESS:
      printf("Con2Prim succeeded.");
      if (set_atmo) {
        printf("Artificial atmosphere has been enforced.");
      }
      if (adjust_cons) {
        printf("Conserved variables have been changed.");
      }
      break;
    case INVALID_DETG:
      printf("Invalid metric determinant: detg = %16.8e ",detg);
      break;
    case NEG_BSQR:
      printf("3-metric not positive, negative B^2: bsq = %16.8e ", bsq);
      break;
    case NANS_IN_CONS:
      printf( "NAN in conserved variables, 
         dens = %16.8e, bsqr = %16.8e, ye = %16.8e", dens, bsqr, ye);
      break;
    case RANGE_RHO:
      printf("Density out of range, dens = %16.8e, rho = %16.8e", dens, rho);
      break;
    case RANGE_EPS:
      printf("Specific energy out of range, eps = %16.8e", eps);
      break;
    case SPEED_LIMIT:
      printf("Speed limit exceeded, v = %16.8e", vel);
      break;
    case B_LIMIT:
      printf("Limit for magnetic field exceeded, b = %16.8e", sqrt(bsqr));
      break;
    case RANGE_YE:
      printf("Electron fraction out of range, Y_e = %16.8e", ye);
      break;
    case ROOT_FAIL_CONV:
      printf("Root finding failed (not converged after %i steps)", iters;
      break;
    case ROOT_FAIL_BRACKET:
      printf("Root finding failed (faulty bracketing)");
      break;
    case PREP_ROOT_FAIL_CONV:
      printf("Preparatory root finding failed (not converged)");
      break;
    case PREP_ROOT_FAIL_BRACKET:
      printf("Preparatory root finding failed (faulty bracketing)");
      break;
    default:
      assert(false);
      printf("Invalid error type. Should never happen. Code is messed up.");
      break;
  }   
} 
  
}; //struct

} //namespace
#endif
