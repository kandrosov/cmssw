def set_gem_alignment_corrections(process,
                                  me11_files=None,
                                  me21_files=None):
    if me11_files is None:
        me11_files = [
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME11_positive_endcap.txt",
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME11_negative_endcap.txt",
        ]
    if me21_files is None:
        me21_files = [
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME21_positive_endcap.txt",
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME21_negative_endcap.txt",
        ]
    process.CSCL1TPLookupTableEP.gemAlignCorrME11Files = me11_files
    process.CSCL1TPLookupTableEP.gemAlignCorrME21Files = me21_files
    return process


def set_gem_alignment_corrections_4d(process,
                                     me11_files=None,
                                     me21_files=None):
    """Use 4D alignment LUTs indexed by (chamber, roll, layer, pad).
    ME11: 36 chambers x 8 rolls x 2 layers x 192 pads = 110592 entries.
    ME21: 18 chambers x 16 rolls x 2 layers x 384 pads = 221184 entries.
    """
    if me11_files is None:
        me11_files = [
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME11_positive_endcap_4d.txt",
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME11_negative_endcap_4d.txt",
        ]
    if me21_files is None:
        me21_files = [
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME21_positive_endcap_4d.txt",
            "L1Trigger/CSCTriggerPrimitives/data/GEMCSC/AlignmentCorrection/GEMCSCLUT_align_corr_es_ME21_negative_endcap_4d.txt",
        ]
    process.CSCL1TPLookupTableEP.gemAlignCorrME11Files = me11_files
    process.CSCL1TPLookupTableEP.gemAlignCorrME21Files = me21_files
    return process
