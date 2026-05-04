#ifndef CondFormats_CSCObjects_CSCL1TPLookupTableME11ILT_h
#define CondFormats_CSCObjects_CSCL1TPLookupTableME11ILT_h

#include "CondFormats/Serialization/interface/Serializable.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <vector>
#include <string_view>
#include <initializer_list>
#include <utility>

#define DECLARE_CSCL1TP_LUT(NAME)    \
private:                             \
  std::vector<unsigned> NAME##_;     \
                                     \
public:                              \
  void set_##NAME(t_lut lut);        \
  unsigned NAME(unsigned idx) const; \
  unsigned NAME##_bit_width() const

#define DEFINE_CSCL1TP_LUT(CLASS, NAME)                                \
  void CLASS::set_##NAME(t_lut lut) { NAME##_ = std::move(lut); }      \
  unsigned CLASS::NAME(unsigned idx) const { return NAME##_.at(idx); } \
  unsigned CLASS::NAME##_bit_width() const { return CSCL1TPLookupTableUtils::get_lut_bit_width(NAME##_); }

struct CSCL1TPLookupTableUtils {
  using t_lut = std::vector<unsigned>;

  static unsigned get_lut_bit_width(const t_lut& lut);
  static unsigned get_common_lut_bit_width(std::initializer_list<unsigned> luts_bit_width,
                                           std::string_view lut_group_name);
  template<unsigned _nChambers, unsigned _nRolls, unsigned _nLayers, unsigned _nPads>
  struct GEMLayout {
    static constexpr unsigned nChambers = _nChambers;
    static constexpr unsigned nRolls = _nRolls;
    static constexpr unsigned nLayers = _nLayers;
    static constexpr unsigned nPads = _nPads;
    static constexpr unsigned size2D = nChambers * nRolls;
    static constexpr unsigned size4D = nChambers * nRolls * nLayers * nPads;

    static void checkIndices(unsigned chamber, unsigned roll, std::optional<unsigned> layer = std::nullopt, std::optional<unsigned> pad = std::nullopt) {
      if (chamber < 1 || chamber > nChambers)
        throw cms::Exception("InvalidLookupTable") << "Chamber number must be in [1," << nChambers << "], but got " << chamber;
      if (roll < 1 || roll > nRolls)
        throw cms::Exception("InvalidLookupTable") << "Roll number must be in [1," << nRolls << "], but got " << roll;
      if (layer.has_value() && (layer.value() < 1 || layer.value() > nLayers))
        throw cms::Exception("InvalidLookupTable") << "Layer number must be in [1," << nLayers << "], but got " << layer.value();
      if (pad.has_value() && pad.value() >= nPads)
        throw cms::Exception("InvalidLookupTable") << "Pad number must be in [0," << nPads - 1 << "], but got " << pad.value();
    }

    static unsigned getIdx4D(unsigned chamber, unsigned roll, unsigned layer, unsigned pad) {
      checkIndices(chamber, roll, layer, pad);
      return nPads * ( nLayers * ( nRolls * (chamber - 1) + (roll - 1) ) + (layer - 1) ) + pad;
    }

    static unsigned getIdx2D(unsigned chamber, unsigned roll) {
      checkIndices(chamber, roll);
      return nRolls * (chamber - 1) + (roll - 1);
    }

    static unsigned getIdx(size_t lut_size, unsigned chamber, unsigned roll, unsigned layer, unsigned pad) {
      if (lut_size == size4D)
        return getIdx4D(chamber, roll, layer, pad);
      else if (lut_size == size2D)
        return getIdx2D(chamber, roll);
      else
        throw cms::Exception("InvalidLookupTable") << "Unsupported LUT size: " << lut_size;
    }

    static void checkLUTSize(size_t size, std::string_view lut_name) {
      if (size != size2D && size != size4D) {
        throw cms::Exception("InvalidLookupTable")
            << lut_name << " size must be either " << size2D << " (2D format) or " << size4D << " (4D format), but got " << size;
      }
    }
  };
};

class CSCL1TPLookupTableME11ILT {
public:
  using t_lut = CSCL1TPLookupTableUtils::t_lut;
  // ME11 layout: chamber in [1,36], roll in [1,8], layer in [1,2], pad in [0,191]
  using GEMLayout = CSCL1TPLookupTableUtils::GEMLayout<36, 8, 2, 192>;

  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME11b_even);
  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME11a_even);
  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME11b_odd);
  DECLARE_CSCL1TP_LUT(GEM_pad_CSC_es_ME11a_odd);

  DECLARE_CSCL1TP_LUT(GEM_roll_CSC_min_wg_ME11_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_CSC_max_wg_ME11_even);
  DECLARE_CSCL1TP_LUT(GEM_roll_CSC_min_wg_ME11_odd);
  DECLARE_CSCL1TP_LUT(GEM_roll_CSC_max_wg_ME11_odd);

  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME11a_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME11a_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME11a_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME11a_odd);

  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME11b_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_2to1_L1_ME11b_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME11b_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_3to1_L1_ME11b_odd);

  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME11a_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME11b_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME11a_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_cosi_corr_L1_ME11b_odd);

  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME11a_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME11b_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME11a_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L1_ME11b_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME11a_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME11b_even);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME11a_odd);
  DECLARE_CSCL1TP_LUT(CSC_slope_corr_L2_ME11b_odd);

  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME11a_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME11a_odd);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME11b_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L1_ME11b_odd);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME11a_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME11a_odd);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME11b_even);
  DECLARE_CSCL1TP_LUT(es_diff_slope_L2_ME11b_odd);

  unsigned es_diff_slope_bit_width() const;

  // GEM alignment corrections in 1/8-strip units.
  // 2D format: chamber x roll.
  // 4D format: chamber x roll x layer x pad.
  // Format is detected automatically from LUT size.
  void set_GEM_align_corr_es_ME11_positive_endcap(std::vector<int> lut);
  void set_GEM_align_corr_es_ME11_negative_endcap(std::vector<int> lut);
  int GEM_align_corr_es_ME11_positive_endcap(unsigned chamber, unsigned roll, unsigned layer, unsigned pad) const;
  int GEM_align_corr_es_ME11_negative_endcap(unsigned chamber, unsigned roll, unsigned layer, unsigned pad) const;

  COND_SERIALIZABLE;

private:
  std::vector<int> GEM_align_corr_es_ME11_positive_endcap_;
  std::vector<int> GEM_align_corr_es_ME11_negative_endcap_;
};

#endif
