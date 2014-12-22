#include "colorramp.hpp"

using namespace std;

//**************************************************************************************************************
unsigned int multi_hue_red_1[27] = {
  255,247,236,  254,232,200,   253,212,158,   253,187,132,   252,141,89,
  239,101,72,   215,48,31,   179,0,0,   127,0,0};

unsigned int multi_hue_red_2[27] = {
  255,255,204,   255,237,160,   254,217,118,   254,178,76,   253,141,60,
  252,78,42,   227,26,28,   189,0,38,   128,0,38};

unsigned int multi_hue_orange[27] = {
  255,255,229,   255,247,188,   254,227,145,   254,196,79,   254,153,41,
  236,112,20,   204,76,2,   153,52,4,   102,37,6};

unsigned int multi_hue_green_1[27] = {
  247,252,253,   229,245,249,   204,236,230,   153,216,201,   102,194,164,
  65,174,118,   35,139,69,   0,109,44,   0,68,27};

unsigned int multi_hue_green_2[27] = {
  255,255,229,   247,252,185,   217,240,163,   173,221,142,   120,198,121,
  65,171,93,   35,132,67,   0,104,55,   0,69,41};

unsigned int multi_hue_bg_1[27] = {
  247,252,240,   224,243,219,   204,235,197,   168,221,181,   123,204,196,
  78,179,211,   43,140,190,   8,104,172,   8,64,129};

unsigned int multi_hue_bg_2[27] = {
  255,247,251,   236,226,240,   208,209,230,   166,189,219,   103,169,207,
  54,144,192,   2,129,138,   1,108,89,   1,70,54};

unsigned int multi_hue_bg_3[27] = {
  255,255,217,   237,248,177,   199,233,180,   127,205,187,   65,182,196,
  29,145,192,   34,94,168,   37,52,148,   8,29,88};

unsigned int multi_hue_blue[27] = {
  255,247,251,   236,231,242,   208,209,230,   166,189,219,   116,169,207,
  54,144,192,   5,112,176,   4,90,141,   2,56,88};

unsigned int multi_hue_purple[27] = {
  247,252,253,   224,236,244,   191,211,230,   158,188,218,   140,150,198,
  140,107,177,   136,65,157,   129,15,124,   77,0,75};

unsigned int multi_hue_pink_1[27] = {
  247,244,249,   231,225,239,   212,185,218,   201,148,199,   223,101,176,
  231,41,138,   206,18,86,   152,0,67,   103,0,31};

unsigned int multi_hue_pink_2[27] = {
  255,247,243,   253,224,221,   252,197,192,   250,159,181,   247,104,161,
  221,52,151,   174,1,126,   122,1,119,   73,0,106};

unsigned int divergent_1[27] = {
  140,81,10,   191,129,45,   223,194,125,   246,232,195,   245,245,245,
  199,234,229,   128,205,193,   53,151,143,   1,102,94};

unsigned int divergent_2[27] = {
  197,27,125,   222,119,174,   241,182,218,   253,224,239,   247,247,247,
  230,245,208,   184,225,134,  127,188,65,   77,146,33};

unsigned int divergent_3[27] = {
  118,42,131,   153,112,171,   194,165,207,   231,212,232,  247,247,247,
  217,240,211,   166,219,160,   90,174,97,   27,120,55};

unsigned int divergent_4[27] = {
  179,88,6,   224,130,20,   253,184,99,   254,224,182,   247,247,247,
  216,218,235,   178,171,210,   128,115,172,   84,39,136};

unsigned int divergent_5[27] = {
  178,24,43,   214,96,77,   244,165,130,   253,219,199,   247,247,247,
  209,229,240,   146,197,222,   67,147,195,   33,102,172};

unsigned int divergent_6[27] = {
  178,24,43,   214,96,77,   244,165,130,   253,219,199,   255,255,255,
  224,224,224,   186,186,186,   135,135,135,   77,77,77};

unsigned int divergent_7[27] = {
  215,48,39,   244,109,67,   253,174,97,   254,224,144,   255,255,191,
  224,243,248,   171,217,233,   116,173,209,   69,117,180};

unsigned int divergent_8[27] = {
  215,48,39,   244,109,67,   253,174,97,   254,224,139,   255,255,191,
  217,239,139,   166,217,106,   102,189,99,   26,152,80};

unsigned int divergent_9[27] = {
  213,62,79,   244,109,67,   253,174,97,   254,224,139,   255,255,191,
  230,245,152,   171,221,164,   102,194,165,   50,136,189};

unsigned int qualitative_1[27] = {
  166,206,227,   31,120,180,   178,223,138,   51,160,44,   251,154,153,
  227,26,28,   253,191,111,   255,127,0,   202,178,214};

unsigned int qualitative_2[27] = {
  251,180,174,   179,205,227,   204,235,197,   222,203,228,   254,217,166,
  255,255,204,   229,216,189,   253,218,236,   242,242,242};

unsigned int qualitative_3[27] = {
  228,26,28,   55,126,184,   77,175,74,   152,78,163,   255,127,0,
  255,255,51,   166,86,40,   247,129,191,   153,153,153};

unsigned int qualitative_4[27] = {
  141,211,199,   255,255,179,   190,186,218,   251,128,114,   128,177,211,
  253,180,98,   179,222,105,   252,205,229,   217,217,217};



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
  current_ramp = static_cast<CRamp>( static_cast<int>(current_ramp) + 1);
  if (current_ramp == CR_END_OF_LIST) current_ramp = (CRamp) 0;
  reset_ramp_array();
  return;
}

void ColorRamp::set_ramp(CRamp ramp_name){
  current_ramp = ramp_name;
  reset_ramp_array();
  return;
}

rgb ColorRamp::get_ramp_color(float norm_value){
  rgb out;

  if (rgb_interpolation){
    cout << "woops: rgb interpolation not yet functional" << endl;
  }
  else{
    out = ramp[(unsigned int)(norm_value*(num_values+1))];
  }

  return out;
}

void ColorRamp::reset_ramp_array(){
  if (ramp != NULL) delete[] ramp;

  unsigned int * rgbvals;
  num_values = 9;
  rgb_max = 255;

  switch (current_ramp){
    case CR_MULTI_HUE_RED_1:
      rgbvals = multi_hue_red_1;
      break;

    case CR_MULTI_HUE_RED_2:
      rgbvals = multi_hue_red_2;
      break;

    case CR_MULTI_HUE_ORANGE:
      rgbvals = multi_hue_orange;
      break;

    case CR_MULTI_HUE_GREEN_1:
      rgbvals = multi_hue_green_1;
      break;

    case CR_MULTI_HUE_GREEN_2:
      rgbvals = multi_hue_green_2;
      break;

    case CR_MULTI_HUE_BG_1:
      rgbvals = multi_hue_bg_1;
      break;

    case CR_MULTI_HUE_BG_2:
      rgbvals = multi_hue_bg_2;
      break;

    case CR_MULTI_HUE_BG_3:
      rgbvals = multi_hue_bg_3;
      break;

    case CR_MULTI_HUE_BLUE:
      rgbvals = multi_hue_blue;
      break;

    case CR_MULTI_HUE_PURPLE:
      rgbvals = multi_hue_purple;
      break;

    case CR_MULTI_HUE_PINK_1:
      rgbvals = multi_hue_pink_1;
      break;

    case CR_MULTI_HUE_PINK_2:
      rgbvals = multi_hue_pink_2;
      break;

    case CR_DIVERGENT_1:
      rgbvals = divergent_1;
      break;

    case CR_DIVERGENT_2:
      rgbvals = divergent_2;
      break;

    case CR_DIVERGENT_3:
      rgbvals = divergent_3;
      break;

    case CR_DIVERGENT_4:
      rgbvals = divergent_4;
      break;

    case CR_DIVERGENT_5:
      rgbvals = divergent_5;
      break;

    case CR_DIVERGENT_6:
      rgbvals = divergent_6;
      break;

    case CR_DIVERGENT_7:
      rgbvals = divergent_7;
      break;

    case CR_DIVERGENT_8:
      rgbvals = divergent_8;
      break;

    case CR_DIVERGENT_9:
      rgbvals = divergent_9;
      break;

    case CR_QUALITATIVE_1:
      rgbvals = qualitative_1;
      break;

    case CR_QUALITATIVE_2:
      rgbvals = qualitative_2;
      break;

    case CR_QUALITATIVE_3:
      rgbvals = qualitative_3;
      break;

    case CR_QUALITATIVE_4:
      rgbvals = qualitative_4;
      break;

    otherwise:
      break;
  }

  ramp = new rgb[num_values];
  for (unsigned int i=0; i<num_values; i++){
    ramp[i].R = float(rgbvals[i*3])/float(rgb_max);
    ramp[i].G = float(rgbvals[i*3+1])/float(rgb_max);
    ramp[i].B = float(rgbvals[i*3+2])/float(rgb_max);
  }

  return;
}
