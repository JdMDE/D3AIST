#ifndef _GAMMATRANSF_H
#define _GAMMATRANSF_H

/**
 @brief Class to implement a gamma transformation of the values of an image, as requested by ITK in the form of a functor
*/
template <class TPixInput, class TPixOutput>
class GammaTransf
{
 public:
  /**
   @brief Default constructor. Do not change.
  */
  GammaTransf() = default;
  /**
   @brief Default destructor. Do not change.
  */
  ~GammaTransf() = default;
  /**
   @brief Procedure to set the gamma value
   \param gamma double The value of gamma to be used
  */
  void SetGammaValue(double gamma) { gamma_value=gamma; };
  /**
   @brief Operator to implement the real trransformation, in this case, normalized_output = (normalized_input)^gamma
  */
  inline TPixOutput operator()(const TPixInput & v) const
  {
    double g;
    g = exp(gamma_value*log(double(v)/mx));
    return TPixOutput(mx*g);
  }
 private:
  double gamma_value;
  double mx=double(itk::NumericTraits<TPixInput>::max());
};

#endif
