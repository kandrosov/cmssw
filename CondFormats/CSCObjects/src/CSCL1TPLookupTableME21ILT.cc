#include "CondFormats/CSCObjects/interface/CSCL1TPLookupTableME21ILT.h"

#define DEFINE_LUT(NAME) DEFINE_CSCL1TP_LUT(CSCL1TPLookupTableME21ILT, NAME)

DEFINE_LUT(GEM_pad_CSC_es_ME21_even);
DEFINE_LUT(GEM_pad_CSC_es_ME21_odd);
DEFINE_LUT(GEM_roll_L1_CSC_min_wg_ME21_even);
DEFINE_LUT(GEM_roll_L1_CSC_max_wg_ME21_even);
DEFINE_LUT(GEM_roll_L1_CSC_min_wg_ME21_odd);
DEFINE_LUT(GEM_roll_L1_CSC_max_wg_ME21_odd);
DEFINE_LUT(GEM_roll_L2_CSC_min_wg_ME21_even);
DEFINE_LUT(GEM_roll_L2_CSC_max_wg_ME21_even);
DEFINE_LUT(GEM_roll_L2_CSC_min_wg_ME21_odd);
DEFINE_LUT(GEM_roll_L2_CSC_max_wg_ME21_odd);
DEFINE_LUT(es_diff_slope_L1_ME21_even);
DEFINE_LUT(es_diff_slope_L1_ME21_odd);
DEFINE_LUT(es_diff_slope_L2_ME21_even);
DEFINE_LUT(es_diff_slope_L2_ME21_odd);
DEFINE_LUT(CSC_slope_cosi_2to1_L1_ME21_even);
DEFINE_LUT(CSC_slope_cosi_2to1_L1_ME21_odd);
DEFINE_LUT(CSC_slope_cosi_3to1_L1_ME21_even);
DEFINE_LUT(CSC_slope_cosi_3to1_L1_ME21_odd);
DEFINE_LUT(CSC_slope_cosi_corr_L1_ME21_even);
DEFINE_LUT(CSC_slope_cosi_corr_L1_ME21_odd);
DEFINE_LUT(CSC_slope_corr_L1_ME21_even);
DEFINE_LUT(CSC_slope_corr_L1_ME21_odd);
DEFINE_LUT(CSC_slope_corr_L2_ME21_even);
DEFINE_LUT(CSC_slope_corr_L2_ME21_odd);

void CSCL1TPLookupTableME21ILT::set_GEM_align_corr_es_ME21_positive_endcap(std::vector<int> lut) {
  GEM_align_corr_es_ME21_positive_endcap_ = std::move(lut);
  GEMLayout::checkLUTSize(GEM_align_corr_es_ME21_positive_endcap_.size(), "GEM_align_corr_es_ME21_positive_endcap");
}
void CSCL1TPLookupTableME21ILT::set_GEM_align_corr_es_ME21_negative_endcap(std::vector<int> lut) {
  GEM_align_corr_es_ME21_negative_endcap_ = std::move(lut);
  GEMLayout::checkLUTSize(GEM_align_corr_es_ME21_negative_endcap_.size(), "GEM_align_corr_es_ME21_negative_endcap");
}

int CSCL1TPLookupTableME21ILT::GEM_align_corr_es_ME21_positive_endcap(unsigned chamber, unsigned roll, unsigned layer, unsigned pad) const {
  const unsigned idx = GEMLayout::getIdx(GEM_align_corr_es_ME21_positive_endcap_.size(), chamber, roll, layer, pad);
  return GEM_align_corr_es_ME21_positive_endcap_.at(idx);
}

int CSCL1TPLookupTableME21ILT::GEM_align_corr_es_ME21_negative_endcap(unsigned chamber, unsigned roll, unsigned layer, unsigned pad) const {
  const unsigned idx = GEMLayout::getIdx(GEM_align_corr_es_ME21_negative_endcap_.size(), chamber, roll, layer, pad);
  return GEM_align_corr_es_ME21_negative_endcap_.at(idx);
}

unsigned CSCL1TPLookupTableME21ILT::es_diff_slope_bit_width() const {
  return CSCL1TPLookupTableUtils::get_common_lut_bit_width(
      {
          es_diff_slope_L1_ME21_even_bit_width(),
          es_diff_slope_L1_ME21_odd_bit_width(),
          es_diff_slope_L2_ME21_even_bit_width(),
          es_diff_slope_L2_ME21_odd_bit_width(),
      },
      "es_diff_slope");
}

#undef DEFINE_LUT
