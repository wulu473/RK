
//! Class to compute 4-th order Runge Kutta update
/** Templated over a type Functor
 *
 */
template<class Functor>
class RK4
{
  public:
    typedef typename Functor::result_type value_type;

    RK4(const Functor& f): m_func(f) {}

    value_type solve(const double t, const value_type& in, const double dt) const
    {
      const value_type k1 = dt*m_func(t        , in        );
      const value_type k2 = dt*m_func(t + dt/2., in + k1/2.);
      const value_type k3 = dt*m_func(t + dt/2., in + k2/2.);
      const value_type k4 = dt*m_func(t + dt   , in + k3   );

      return in + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
    }

  private:
    Functor m_func;
};

