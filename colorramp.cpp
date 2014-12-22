#include "colorramp.hpp"

using namespace std;

ColorRamp::ColorRamp(){
  ramp = NULL;
  rgb_max = 255;
  num_values = 9;
  rgb_interpolation = false;
  current_ramp = CR_DIVERGENT_8;
  reset_ramp_array();
}

ColorRamp::~ColorRamp(){
  if (ramp != NULL) delete[] ramp;
}

void ColorRamp::cycle_ramp(){
  current_ramp++;
  if (current_ramp == CR_END_OF_LIST) current_ramp = 0;
  reset_ramp_array();
  return;
}

void ColorRamp::set_ramp(CRamp ramp_name){
  current_ramp = ramp_name;
  reset_ramp_array();
  return;
}

rgb ColorRamp::get_ramp_color(float norm_value){
  if (rgb_interpolation){

  }
  else return ramp[unsigned int(norm_value*(num_values+1))];
}

void ColorRamp::reset_ramp_array(){
  if (ramp != NULL) delete[] ramp;

  unsigned int * rgbvals;
  num_values = 9;
  rgb_max = 255;

  switch current_ramp{
    case CR_MULTI_HUE_RED_1:
      rgbvals = multi_hue_red_1;
      return;

    case CR_MULTI_HUE_RED_2:
      rgbvals = multi_hue_red_2;
      return;

    case CR_MULTI_HUE_ORANGE:
      rgbvals = multi_hue_orange;
      return;

    case CR_MULTI_HUE_GREEN_1:
      rgbvals = multi_hue_green_1;
      return;

    case CR_MULTI_HUE_GREEN_2:
      rgbvals = multi_hue_green_2;
      return;

    case CR_MULTI_HUE_BG_1:
      rgbvals = multi_hue_bg_1;
      return;

    case CR_MULTI_HUE_BG_2:
      rgbvals = multi_hue_bg_2;
      return;

    case CR_MULTI_HUE_BG_3:
      rgbvals = multi_hue_bg_3;
      return;

    case CR_MULTI_HUE_BLUE:
      rgbvals = multi_hue_blue;
      return;

    case CR_MULTI_HUE_PURPLE:
      rgbvals = multi_hue_purple;
      return;

    case CR_MULTI_HUE_PINK_1:
      rgbvals = multi_hue_pink_1;
      return;

    case CR_MULTI_HUE_PINK_2:
      rgbvals = multi_hue_pink_2;
      return;

    case CR_DIVERGENT_1:
      rgbvals = divergent_1;
      return;

    case CR_DIVERGENT_2:
      rgbvals = divergent_2;
      return;

    case CR_DIVERGENT_3:
      rgbvals = divergent_3;
      return;

    case CR_DIVERGENT_4:
      rgbvals = divergent_4;
      return;

    case CR_DIVERGENT_5:
      rgbvals = divergent_5;
      return;

    case CR_DIVERGENT_6:
      rgbvals = divergent_6;
      return;

    case CR_DIVERGENT_7:
      rgbvals = divergent_7;
      return;

    case CR_DIVERGENT_8:
      rgbvals = divergent_8;
      return;

    case CR_DIVERGENT_9:
      rgbvals = divergent_9;
      return;

    case CR_QUALITATIVE_1:
      rgbvals = qualitative_1;
      return;

    case CR_QUALITATIVE_2:
      rgbvals = qualitative_2;
      return;

    case CR_QUALITATIVE_3:
      rgbvals = qualitative_3;
      return;

    case CR_QUALITATIVE_4:
      rgbvals = qualitative_4;
      return;

    otherwise:

  }

  ramp = new rgb[num_values];
  for (unsigned int i=0; i<num_values; i++){
    ramp[i].R = float(rgbvals[i*3])/float(rgb_max);
    ramp[i].G = float(rgbvals[i*3+1])/float(rgb_max);
    ramp[i].B = float(rgbvals[i*3+2])/float(rgb_max);
  }

  return;
}
