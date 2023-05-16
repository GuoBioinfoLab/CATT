const ref_prefix = PROGRAM_FILE[1:end-7] * "resource"
const bwa_path = "bwa"
const samtools_path = "/Users/kroaity/Documents/bin/samtools"

const refg = Dict(

    "hs" => Dict(

        "TRA" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRA/TRAV-gai3-DNA.fa",
                "jregion" => "$ref_prefix/TR/hs/TRA/TRAJ-gai2-DNA.fa",
                "cmotif" => "[EHILMSTV]{1}Y[FILY]{1}C[AGILV]{1}",
                "fmotif" => "(LA|YI|FI|II|LY|LM|PT|TI|LV|ST|VT|LT|LI|LQ|MR|VI|FV|FQ|LF|LL|FE|FT|LS|LN|FY)F((ARG)|(G[A-Z]{1}G))",
                "coffset" => 3,
                "foffset" => 0,
                "innerC" => "place_holder",
                "innerF" => "place_holder",
            )
        ),

        "TRB" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRB/TRBV-gai6-DNA.fa",
                "jregion" => "$ref_prefix/TR/hs/TRB/TRBJ-gai3-DNA.fa",
                "cmotif" => "(LR|YF|YI|YL|YQ|YR)C(A|S|T|V|G|R|P|D)",
                "fmotif" => "[FTYH]{1}FG[ADENPQS]{1}G",
                "coffset" => 2,
                "foffset" => 1,
                "innerC" => "(CAS|CSA|CAW|CAT|CSV|CAI|CAR)",
                "innerF" => "place_holder",
                "Dregion" => [
                    ("TRBD1*01", "gggacagggggc"),
                    ("TRBD2*01", "gggactagcggggggg"),
                    ("TRBD2*02", "gggactagcgggaggg")
                ]
            )
        ),

        "TRG" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRGV.fa",
                "jregion" => "$ref_prefix/TR/hs/TRGJ.fa",
                "cmotif" => "(YY|YH)C(A|T)",
                "fmotif" => "([LV]{1}FG[SP]{1}G)|(IFAEG)|(TFAKG)",
                "coffset" => 2,
                "foffset" => 1,
                "innerC" => "(CAT|CAW|CTT|CAL|CAC|CAA)",
                "innerF" => "(KLF)|(KVF)|(KIF)|(KTF)",
            )
        ),

        "TRD" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRDV.fa",
                "jregion" => "$ref_prefix/TR/hs/TRDJ.fa",
                "cmotif"  => "(YF|YL|YY)C(A|S)",
                "fmotif"  => "[I|F]FG[K|T]G",
                "coffset" => 2,
                "foffset" => 1,
                "innerC"  => "place_holder",
                "innerF"  => "place_holder",
                "Dregion" => [
                    ("TRDD1*01", "gaaatagt"),
                    ("TRDD2*01", "ccttcctac"),
                    ("TRDD3*01", "actgggggatacg")
                ]
            )
        ),

        "IGH" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/IG/hs/IGHV-gai2-DNA.fa",
		        "cmotif" => "(LYY|TAV|YVA|TYY|MYY|VYD|AFN|AYY|EYY|THY|VYY|LYH|VDY)C",
                "jregion" => "$ref_prefix/IG/hs/IGHJ.fa",
                "fmotif"  => "[VLHYPIS]{1}WG[QR]{1}G",
                "coffset" => 3,
                "foffset" => 1,
                "innerC" => "place_holder",
                "innerF" => "place_holder",
            ),
        ),
    ),

    "ms" => Dict(
        "TRB" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/ms/TRBV-CDR3.fa",
                "jregion" => "$ref_prefix/TR/ms/TRBJ-CDR3.fa",
                "cmotif" => "[AGPST]{1}C(LHV|FYI|LYF|LFV|FYV|LYL|MYL|LYT|TCY|LYV|FYL|FYT|LYM|LYI|FYM)",
                "fmotif" => "(SF|FF|YF|HF|TF|DH|LF)G[ADMREKLSPH]{1}G",
                "innerC" => "place_holder",
                "innerF" => "place_holder",
                "coffset" => 0,
                "foffset" => 1,
            ),

            # "CDR2" => Dict(),
        ),
        # "TRA" => Dict(),
        # "IGH" => Dict(),
        # "IGK" => Dict()
    ),

    "pig" => Dict(
        "TRB" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/TR/pig/TRBV-gai1-DNA.fa",
                "cmotif"  => "(FF|FL|LF|YF|YL|YV)C",
                "coffset" => 2,
                "innerC"  => "place_holder",

                "jregion" => "$ref_prefix/TR/pig/TRBJ-gai1-DNA.fa",
                "innerF"  => "(LH|LI|LT|LY|QH|QI|QY|RY|VF|YN)F",
                "fmotif"  => "[HILTYFM]{1}FG[ADEGLPS]{1}G",
                "foffset" => 1,

                "Dregion" => [
                    ("TRBD1*01", "gggacagggggggc"),
                    ("TRBD2*01", "ggagctatgggggggg"),
                    ("TRBD3*01", "ggagctatggggggggg")
                ]
            )
        )
    )
)
