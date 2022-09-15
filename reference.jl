const ref_prefix = PROGRAM_FILE[1:end-12] * "resource"
const bwa_path = "bwa"
const samtools_path = "/home/chensy/anaconda3/bin/samtools"

const refg = Dict(

    "hs" => Dict(
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
            ),
            "CDR1" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRB/TRBV-CDR1-front.fa",
                "jregion" => "$ref_prefix/TR/hs/TRB/TRBV-CDR1-back.fa",
                "cmotif" => "(HQT|KPI|SQT|YPI|VPI|QPI|API|TPK|DSI|SPI|NPI|EPI|RSL|QVD|TPE|DPS|EQH|TQD|DPI|TVE|AQD|SSQ|SPK|SPM|PQN|VQD|SQN|EQN|SPR)C",
                "fmotif" => "(IF|MY|MF|VF|VH|FY|LY|VY|LL|MS|VC|LS|LF|VS)|W",
                "coffset" => 3,
                "foffset" => 2,
                "innerC" => "place_holder",
                "innerF" => "place_holder",
            ),
            "CDR2" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRB/TRBV-CDR2-front.fa",
                "jregion" => "$ref_prefix/TR/hs/TRB/TRBV-CDR2-back.fa",
                "cmotif" => "LLW|LTY|LIR|LVS|LNY|IYY|IYF|LIS|LVY|LIH|LFE|LLH|LVQ|LFY|MFV|FIY|LIY|VYY|LIQ|IHY|IFQ|MAT|LAY|IFE|MAA|MVY|LSY|MFS|FIS|LFH|IAT",
                "fmotif" => "NNE|INK|IDE|AVD|QIS|RAK|RNR|FQK|IIN|LEK|KEK|SDK|TAK|NGR|LDD|TGK|QDE|SIN|QDK|RDK|RQR|TDK|RNK|ADD|RGR|TEN|IEK|IQK|FDE|IDD|PDK|LDK|RSE|SEK|TDQ|VEN|IDN|TYE|TEK|TNK|LQE|VDD",
                "innerC" => "place_holder",
                "innerF" => "place_holder",
                "coffset" => 0,
                "foffset" => 1,
            )
        ),

        "TRA" => Dict(

            "CDR1" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRA/TRAV-CDR1-front.fa",
                "jregion" => "$ref_prefix/TR/hs/TRA/TRAJ-CDR1-front.fa",
                "cmotif" => "C(TSS|SYK|TYD|TYS|NHS|SFT|NSS|NFS|SSS|TYT|VYE|AYE|SYE|SYT|SYS|TYQ|SCT|SFP|NYS|AYS|DYT|SHN|NYT)",
                "fmotif" => "(IQ|LH|FF|WL|FM|FP|LQ|VH|IH|LR|VY|LY|LS|FL|LL|FI|FH|LF|IT)|W",
                "coffset" => -1,
                "foffset" => 2,
                "innerC" => "place_holder",
                "innerF" => "place_holder",
            ),
            "CDR2" => Dict(
                "vregion" => "$ref_prefix/TR/hs/TRA/TRAV-CDR2-front.fa",
                "jregion" => "$ref_prefix/TR/hs/TRA/TRAV-CDR2-back.fa",
                "cmotif" => "LTY|LIR|LTI|LMM|LIL|LIA|IIH|LIS|LFM|LVT|LQR|IID|LTL|VIH|LFY|LLL|LIY|LIQ|VIR|IIQ|IMF|LMS|LSY|LMI|LLR|IMS|LFT|LMY|HLK|LLS|LVK|LLK|LFV",
                "fmotif" => "ETN|KKQ|KQE|KGH|KSH|KHS|QTS|KEK|KRH|EIS|GRN|KKD|ESI|KQD|KQN|KKS|MRR|KGI|KED|QGI|ATE|KSN|EDG|EEK|KRK|NSK|TSN|KQK|GSN|KKH|EET|TDS|KKK|KKE|VTN|QGD|NED|KGS|KKL|SQQ|VNN|RQG",
                "innerC" => "place_holder",
                "innerF" => "place_holder",
                "coffset" => 0,
                "foffset" => 1,
            ),
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

        "IGH" => Dict(
            "CDR3" => Dict(
                "vregion" => "$ref_prefix/IG/hs/IGHV-gai2-DNA.fa",
		        "cmotif" => "(LYY|TAV|YVA|TYY|MYY|VYD|AFN|AYY|EYY|THY|VYY|LYH|VDY)C",
                "jregion" => "$ref_prefix/IG/hs/IGHJ.fa",
                "fmotif"  => "[VLHYPIS]{1}WG[QR]{1}G",
                "coffset" => 3,
                "foffset" => 1,
                "innerC" => "ASDASDASDASDADS",
                "innerF" => "ASDASDASDASDASD",
            ),
        ),

        # "IGK" => Dict(

        # )
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
        ),
        # "TRA" => Dict(
        #     "CDR3" => Dict(
        #         "vregion" => "ASD",
        #     )
        # )
    )

)
