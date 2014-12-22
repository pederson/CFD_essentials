#ifndef COLORRAMP_H
#define COLORRAMP_H

#include <stdlib.h>
#include <iostream>

enum CRamp{
      CR_MULTI_HUE_RED_1 = 0, CR_MULTI_HUE_RED_2, CR_MULTI_HUE_ORANGE, CR_MULTI_HUE_GREEN_1, CR_MULTI_HUE_GREEN_2,
      CR_MULTI_HUE_BG_1, CR_MULTI_HUE_BG_2, CR_MULTI_HUE_BG_3, CR_MULTI_HUE_BLUE, CR_MULTI_HUE_PURPLE,
      CR_MULTI_HUE_PINK_1, CR_MULTI_HUE_PINK_2,
      CR_DIVERGENT_1, CR_DIVERGENT_2, CR_DIVERGENT_3, CR_DIVERGENT_4, CR_DIVERGENT_5, CR_DIVERGENT_6,
      CR_DIVERGENT_7, CR_DIVERGENT_8, CR_DIVERGENT_9,
      CR_QUALITATIVE_1, CR_QUALITATIVE_2, CR_QUALITATIVE_3, CR_QUALITATIVE_4,
      CR_END_OF_LIST};


struct rgb{
  float R;
  float G;
  float B;
};

class ColorRamp{
public:

  ColorRamp();
  ~ColorRamp();

  void cycle_ramp();
  void set_ramp(CRamp ramp_name);
  rgb get_ramp_color(float norm_value);

private:

  rgb * ramp;
  CRamp current_ramp;
  unsigned int rgb_max;
  unsigned int num_values;
  bool rgb_interpolation;

  void reset_ramp_array();

};


#endif