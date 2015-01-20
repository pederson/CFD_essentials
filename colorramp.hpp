#ifndef COLORRAMP_H
#define COLORRAMP_H

#include <stdlib.h>
#include <iostream>

enum class CRamp : char {
      MULTI_HUE_RED_1, MULTI_HUE_RED_2, MULTI_HUE_ORANGE, MULTI_HUE_GREEN_1, MULTI_HUE_GREEN_2,
      MULTI_HUE_BG_1, MULTI_HUE_BG_2, MULTI_HUE_BG_3, MULTI_HUE_BLUE, MULTI_HUE_PURPLE,
      MULTI_HUE_PINK_1, MULTI_HUE_PINK_2,
      DIVERGENT_1, DIVERGENT_2, DIVERGENT_3, DIVERGENT_4, DIVERGENT_5, DIVERGENT_6,
      DIVERGENT_7, DIVERGENT_8, DIVERGENT_9,
      QUALITATIVE_1, QUALITATIVE_2, QUALITATIVE_3, QUALITATIVE_4,
      END_OF_LIST};


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
  rgb get_ramp_color(float norm_value) const;

private:

  rgb * ramp;
  CRamp current_ramp;
  unsigned int rgb_max;
  unsigned int num_values;
  bool rgb_interpolation;

  void reset_ramp_array();

};

class Color{
public:

  static rgb Red();
  static rgb Blue();
  static rgb Green();
  static rgb Orange();
  static rgb Yellow();
  static rgb Violet();
  static rgb Pink();
  static rgb White();
  static rgb Black();
  static rgb Brown();

};


#endif