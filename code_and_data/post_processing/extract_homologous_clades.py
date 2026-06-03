import os

from Bio import SeqIO
from ete3 import Tree


def filter_fasta_by_ids(fasta_file, id_list, output_file):
    """Read a FASTA file, keep only sequences whose ID is in ``id_list``, and write them out."""
    sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    with open(output_file, 'w') as fasta_output:
        for seq_id in id_list:
            if seq_id in sequence_dict:
                fasta_output.write(f">{seq_id}\n{sequence_dict[seq_id]}\n")
            else:
                print(f"ID {seq_id} not found in the FASTA file.")


def find_leaf_nodes_within_boundaries(boundaries, newick_file):
    """Collect leaves under each boundary's MRCA.

    Boundary formats:
      [left_id, right_id]      -> include all leaves under the clade
      [left_id, right_id, -1]  -> exclude all leaves under the clade
    """
    tree = Tree(newick_file, format=1)

    include_list = []
    exclude_set = set()

    for boundary in boundaries:
        if len(boundary) == 3 and boundary[2] == -1:
            left_id, right_id, _ = boundary
            mode = "exclude"
        elif len(boundary) == 2:
            left_id, right_id = boundary
            mode = "include"
        else:
            raise ValueError(
                f"Invalid boundary spec {boundary!r}: expected "
                "[left, right] or [left, right, -1]."
            )

        current_node = tree.get_common_ancestor(left_id, right_id)
        current_subtree_leaf_names = current_node.get_leaf_names()

        if mode == "include":
            include_list += current_subtree_leaf_names
        else:
            exclude_set.update(current_subtree_leaf_names)

    # De-duplicate while preserving insertion order, dropping anything in exclude_set.
    seen = set()
    unique_list = []
    excluded_count = 0
    for name in include_list:
        if name in exclude_set:
            excluded_count += 1
            continue
        if name not in seen:
            unique_list.append(name)
            seen.add(name)

    if exclude_set:
        print(f"Excluded {excluded_count} leaves under {len(exclude_set)} excluded-clade leaf names.")
    print("Total number of identified orthologs:", len(unique_list))
    return unique_list


def prune_tree_by_ids(tree_file, keep_ids, output_tree_file):
    """Read a Newick tree, prune to ``keep_ids``, and write the resulting subtree."""
    tree = Tree(tree_file, format=1)

    existing_ids = set(tree.get_leaf_names())
    keep_present = [leaf_id for leaf_id in keep_ids if leaf_id in existing_ids]
    missing = len(keep_ids) - len(keep_present)
    if missing:
        print(f"{missing} ortholog IDs were not found in {os.path.basename(tree_file)}")
    if not keep_present:
        print(f"No ortholog IDs found in {tree_file}; skipping tree output.")
        return False

    tree.prune(keep_present, preserve_branch_length=True)
    tree.write(outfile=output_tree_file)
    return True


def derive_hmm_tree_paths(gene_name, tree_dir):
    """Return the two candidate HMM-ordered tree filenames for a gene (regions / no-regions)."""
    return [
        os.path.join(tree_dir, f"{gene_name}_hmm_E1000_db_hmmregions_FAMSA_gt0.1_hmmordered.tree"),
        os.path.join(tree_dir, f"{gene_name}_hmm_E1000_db_FAMSA_gt0.1_hmmordered.tree"),
    ]


def derive_m8_ortholog_tree_path(gene_name, tree_file, tree_dir):
    """Build the output ortholog-tree path for an m8-ordered input tree."""
    suffix = os.path.basename(tree_file)
    if "_db_" in suffix:
        suffix = suffix.split("_db_", 1)[1]
        return os.path.join(tree_dir, f"{gene_name}_db_{suffix.replace('.tree', '_orthologs.tree')}")
    return os.path.join(tree_dir, f"{gene_name}_{suffix.replace('.tree', '_orthologs.tree')}")


# -----------------------------------------------------------------------------
# USER CONFIGURATION
# Replace each path below with the corresponding location on your machine.
# `HMM_TREE_DIR`     : folder containing the HMM- and m8-ordered Newick trees.
# `MSA_DIR`          : folder containing the tree-ordered FASTA alignments.
# `ORTHOLOG_OUT_DIR` : folder where the per-gene ortholog FASTAs will be written.
# -----------------------------------------------------------------------------
HMM_TREE_DIR     = r"/path/to/hmmorder_trees"
MSA_DIR          = r"/path/to/treeorder_msa"
ORTHOLOG_OUT_DIR = r"/path/to/ortholog_lists"

#%% Gene boundaries
gene_boundaries = {
    "CsrA": {"boundary": [["NZ_CP018846.1_3238", "JAGBIB010000205.1_6"]], "hmmregion": 1}, # Done
    "FlaG": {"boundary": [["CAKRAK010000006.1_63","NZ_LR884459.1_835"]], "hmmregion": 0}, # Done
    "FlgA": {"boundary": [["NZ_FMWO01000048.1_100", "PLHP01000069.1_11"]], "hmmregion": 0}, # Done
    "FlgB": {"boundary": [["NZ_SMCO01000001.1_696", "JAABRT010000039.1_14"]], "hmmregion": 0}, # Done
    "FlgC": {"boundary": [["NZ_QJKC01000011.1_38", "WQVC01000133.1_2"]], "hmmregion": 0}, # Done
    "FlgD": {"boundary": [["NZ_JRUD01000022.1_3", "JAGDWZ010000024.1_35"]], "hmmregion": 1}, # Done
    "FlgE": {"boundary": [["NZ_LXUF01000016.1_129", "DZAY01000028.1_17"]], "hmmregion": 0}, # Done
    "FlgF": {"boundary": [["NZ_SNWX01000003.1_162", "SHAP01000004.1_89"]], "hmmregion": 0}, # Done
    "FlgG": {"boundary": [["JAAZZB010000056.1_7","JACRGQ010000032.1_11"]],"hmmregion": 0}, # Done
    "FlgH": {"boundary": [["NZ_KE384096.1_251", "LMZT01000109.1_25"],["NZ_CP029161.1_246","NZ_CP034852.1_251"]], "hmmregion": 0}, # Done
    "FlgI": {"boundary": [["JAGXAQ010000001.1_246", "JAHBXO010000001.1_195"]], "hmmregion": 0}, # Done
    "FlgJ": {"boundary": [["NZ_CM001773.1_1530","VGJY01000006.1_37"]] , "hmmregion": 0}, # Done
    "FlgK": {"boundary": [["NZ_LT984884.1_1377", "NZ_LFLZ01000002.1_179"]], "hmmregion": 0}, # Done
    "FlgL": {"boundary": [["JAENYK010000004.1_230","JABSOS010000006.1_147"]], "hmmregion": 0}, # Done
    "FlgM": {"boundary": [["AP019373.1_619", "JAFMML010000021.1_46"]], "hmmregion": 1}, # Done
    "FlgO": {"boundary": [["JABSRV010000072.1_8", "CAMDPZ010000038.1_12"]], "hmmregion": 0}, # Done
    "FlgP": {"boundary": [["NZ_LFLS01000069.1_1", "CP071184.1_532"]], "hmmregion": 1}, # Done
    "FlgT": {"boundary": [["NZ_WOBR01000004.1_5", "JAFGUH010000219.1_18"],["JAGOXX010000126.1_1","DRBZ01000012.1_5"]], "hmmregion": 1},
    "FlgQ": {"boundary": [["NZ_CP043427.1_1608", "DQOT01000262.1_1"]], "hmmregion": 0}, # Done
    "FlhA": {"boundary": [["VFJK01000130.1_22", "JAAZXR010000010.1_48"]], "hmmregion": 0}, # Done
    "FlhB": {"boundary": [["JADGAT010000004.1_80", "PWFS01000245.1_3"]], "hmmregion": 0}, # Done
    "FlhC": {"boundary": [["NZ_CADCXY010000003.1_57", "JABTUS010000012.1_13"]], "hmmregion": 0}, # Done
    "FlhD": {"boundary": [["NZ_FCOP01000013.1_300", "DVOQ01000054.1_2"], ["DBVL01000089.1_7", "NZ_MCAS01000001.1_498"], 
                          ["NZ_CP026111.1_3023", "NZ_QKOB01000016.1_19"]], "hmmregion": 1}, # Done
    "FlhE": {"boundary": [["NZ_PUJU01000006.1_43", "NZ_LYRP01000001.1_468"], 
                          ["NZ_AQWM01000013.1_60", "NZ_LFLS01000010.1_141"]], "hmmregion": 0}, #Done
    "FlhF": {"boundary": [["JAFIFN010000118.1_1", "DCTH01000103.1_6"]], "hmmregion": 0}, # Done
    "FlhG": {"boundary": [["NZ_CP074572.1_695", "NZ_CP036289.1_1414"]], "hmmregion": 1}, # Done
    "FliA": {"boundary": [[ "NZ_PVNH01000004.1_217","JACQDB010000019.1_3"]], "hmmregion": 0}, # Done
    "FliC": {"boundary": [["JABTUQ010000007.1_105", "NZ_NPEV01000072.1_21"]], "hmmregion": 1}, # Done
    "FliD": {"boundary": [["VGJQ01000016.1_9", "NZ_RQYL01000006.1_16"],["JAKAVV010000005.1_147","JAJFBP010000037.1_1"],
                          ["CP070841.1_1245","NZ_FOJI01000003.1_297"]], "hmmregion": 0}, # Done
    "FliE": {"boundary": [["NZ_JACHWY010000001.1_598","NZ_NCXP01000045.1_13"]], "hmmregion": 0}, # Done
    "FliF": {"boundary": [["KI928935.1_128", "NZ_NFMA01000022.1_40"], ["JALNZK010000024.1_7", "NZ_CP020612.1_2605"]], "hmmregion": 0}, # Done
    "FliG": {"boundary": [["QMQA01000039.1_3","WGDQ01000488.1_4"]], "hmmregion": 0}, # Done
    "FliH": {"boundary": [["NZ_KB905155.1_716", "NZ_NNRM01000020.1_27"]], "hmmregion": 1}, # Done
    "FliI": {"boundary": [["NZ_CP014839.1_3861","NZ_JAMXSR010000002.1_3082" ]], "hmmregion": 0}, # Done 
    "FliJ": {"boundary": [["NZ_CP017705.1_377", "JADLDL010000282.1_8"]], "hmmregion": 1}, # Done
    "FliK": {"boundary": [["NZ_JACDRR010000011.1_72","SCUE01000012.1_199"],
                          ["JAAEQF010000458.1_9","NZ_CCDG010000024.1_25"]],"hmmregion": 1}, # Done
    "FliL": {"boundary": [["JAGOGO010000029.1_5","JALHOL010000001.1_89"]], "hmmregion": 1}, # Done
    "FliM": {"boundary": [["WJPL01000003.1_718","NZ_RAVW01000244.1_8"]], "hmmregion": 0}, # Done
    "FliN": {"boundary": [["DDHO01000024.1_3", "NZ_JANUCP010000003.1_99"]], "hmmregion": 1}, # Done
    "FliO": {"boundary": [["DPKF01000115.1_9", "JAILEI010000204.1_40"]], "hmmregion": 1}, # Done
    "FliP": {"boundary": [["PMZA01000392.1_5", "JACAED010000054.1_36"]], "hmmregion": 0}, # Done
    "FliQ": {"boundary": [["NZ_CAJZAH010000004.1_271", "NZ_CP008852.1_2500"]], "hmmregion": 0}, # Done
    "FliR": {"boundary": [["JAJFTF010000051.1_51", "CP022114.1_124"]], "hmmregion": 0}, # Done
    "FliS": {"boundary": [["CADEEX010000001.1_29", "WRFW01000116.1_5"]], "hmmregion": 1}, # Done 
    "FliW": {"boundary": [["JAGPBS010000023.1_7", "JAIUNT010000009.1_15"], ["CAKVCL010000012.1_9", "NZ_CP035928.1_1262"]], "hmmregion": 0}, # Done 
    "FliZ": {"boundary": [["NZ_CABHYG010000009.1_70", "UGTZ01000002.1_29"]], "hmmregion": 0}, # Done 
    "MotA": {"boundary": [["NZ_LROS01000079.1_43","WQZS01000003.1_66"], ["NZ_JQAR01000006.1_92", "JALZUD010000138.1_7"]], "hmmregion": 0},  # Done
    "MotB": {"boundary": [["NZ_QRAP01000001.1_885", "JABDKO010000122.1_6"]], "hmmregion": 1}, # Done
    "SwrD": {"boundary": [["NZ_LWAF01000001.1_42", "NZ_WTLI01000012.1_102"]], "hmmregion": 0}, # Done
    "FlgN": {"boundary": [["NC_014121.1_3232","JAAOMR010000134.1_46"],["PWLT01000094.1_3","CAJJBF010000012.1_15"],["JAAEMS010000170.1_50","JAHDQN010000020.1_199"]],
             "hmmregion": 1}, # Done
    "FliT": {"boundary": [["JAMDYG010000001.1_1082","NZ_LHOX01000021.1_174"]],"hmmregion":1}, # Done
    "FlbT": {"boundary": [["NZ_VHLH01000025.1_23","NZ_NHSK01000154.1_13"]], "hmmregion": 1}, # Done
    "FlaF": {"boundary": [["NZ_RJRZ01000020.1_48","QHBU01000122.1_3"]],"hmmregion": 0}, # Done
    "Putative": {"boundary": [["JALOAD010000056.1_9","DIVV01000130.1_6"]],"hmmregion":0}, #Done
    "PilZ": {"boundary": [["NZ_MSLT01000023.1_561","JAAYCS010000158.1_25"]],"hmmregion":0}, #Done
    "MotE": {"boundary": [["NZ_AP024401.1_284","NZ_FRAU01000010.1_67"],["JAEYOK010000016.1_113","JAILQQ010000206.1_11"],["CAMAQH010000001.1_43","JACPIN010000018.1_13"]],"hmmregion":0}, #Done
    "MotY": {"boundary": [["NZ_KB822603.1_18","CAIRAY010000119.1_4"]],"hmmregion":0}, #Done
    "SwrB": {"boundary": [["JAFLRU010000043.1_11","NZ_CP048020.1_146"]],"hmmregion":0,
             "tree": fr"{HMM_TREE_DIR}\SwrB_db_FAMSA_gt0.1_m8ordered_rerooted.tree",
             "fasta": fr"{MSA_DIR}\SwrB_db_FAMSA_gt0.1_treeordered.fasta"}, # Done
    "PflA": {"boundary": [["NZ_LN831025.1_1648","MAAG01000128.1_20"]],"hmmregion":0,
             "tree": fr"{HMM_TREE_DIR}\PflA_db_FAMSA_gt0.1_m8ordered.tree",
             "fasta": fr"{MSA_DIR}\PflA_db_FAMSA_gt0.1_treeordered.fasta"}, # Done
    "PflB": {"boundary": [["NZ_LN831025.1_380","JAIOSF010000021.1_111"]],"hmmregion":0,
             "tree": fr"{HMM_TREE_DIR}\PflB_db_FAMSA_gt0.1_m8ordered.tree",
             "fasta": fr"{MSA_DIR}\PflB_db_FAMSA_gt0.1_treeordered.fasta"}, # Done
    "MotC": {"boundary": [["NZ_JACHIK010000011.1_65","JAHDFM010000008.1_6"]],"hmmregion":0,
             "tree": fr"{HMM_TREE_DIR}\MotC_db_FAMSA_gt0.1_m8ordered.tree",
             "fasta": fr"{MSA_DIR}\MotC_db_FAMSA_gt0.1_treeordered.fasta"}, # Done
    "FliB": {"boundary": [["CAKQYU010000002.1_226","CAKQKE010000010.1_26"]],"hmmregion":0}, #Done 
    "MotX": {"boundary": [["AP025472.1_2530","CAAGDP010000202.1_7"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\MotX_VIBPA_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\MotX_VIBPA_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "MotK": {"boundary": [["NC_007493.2_2953","NZ_JABFCX010000003.1_989"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\MotK_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\MotK_db_FAMSA_gt0.1_treeordered.fasta"
                }, #Done
    "FlrC": {"boundary": [["NZ_BBIU01000036.1_34","QKED01000117.1_6"]],"hmmregion":0}, #Done
    "FlrA": {"boundary": [["NZ_CP018616.1_1940","JAABTJ010000001.1_62"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlrA_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlrA_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "SwrA": {"boundary": [["NZ_LSBB01000001.1_575","NZ_JAKZKS010000002.1_50"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\SwrA_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\SwrA_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "FlcD": {"boundary": [["NC_001318.1_230","DNNS01000309.1_8"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlcD_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlcD_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "FlcA": {"boundary": [["NC_001318.1_316","NZ_AP024401.1_322"],["QNBN01000028.1_5","JAFGXN010000061.1_1"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlcA_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlcA_db_FAMSA_gt0.1_treeordered.fasta"
                }, #A secondary clade that is also annoated as FlcA but it is not very similar sequence-wise
    "FlcB": {"boundary": [["NC_001318.1_56","JACRPF010000011.1_139"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlcB_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlcB_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "FlcC": {"boundary": [["NZ_JACHFC010000001.1_218","JAKSCK010000028.1_18"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlcC_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlcC_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "FlaY": {"boundary": [["NC_014375.1_2630","JAJZGY010000023.1_20"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FlaY_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FlaY_db_FAMSA_gt0.1_treeordered.fasta"
                },
    "FapA": {"boundary": [["NZ_CP069213.1_501","JAFIGE010000066.1_15"]],"hmmregion":0}, # second clade is MinC
    "DUF1217": {"boundary": [["NZ_MDET01000003.1_55","JACQAI010000194.1_6"]],"hmmregion":0},
    "DUF6470": {"boundary": [["NZ_CABKRX010000062.1_30","CP060226.1_2110"],["JAAYOS010000013.1_8","JAHHUB010000032.1_2"]], "hmmregion":0},
    "DUF327": {"boundary": [["CAJUQB010000006.1_41","JAJXUE010000021.1_9"]], "hmmregion":0},
    "Transglycosylase":{"boundary":[["UCEZ01000055.1_32","DFNM01000065.1_5"]],"hmmregion":0},
    "FljA" : {"boundary":[["NC_003197.2_2734","JAEWTL010000044.1_18"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\FljA_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\FljA_db_FAMSA_gt0.1_treeordered.fasta"},
    "YdiV" : {"boundary":[["NZ_CP040443.1_4484","NZ_JAGGMQ010000001.1_2054"]],"hmmregion":0,
                "tree": fr"{HMM_TREE_DIR}\YdiV_db_FAMSA_gt0.1_m8ordered.tree",
                "fasta": fr"{MSA_DIR}\YdiV_db_FAMSA_gt0.1_treeordered.fasta"},
    "YvyF" : {"boundary":[["NZ_CP026362.1_4267","JAAXZX010000068.1_1"]],"hmmregion":0},
    "FlgR" : {"boundary":[["NC_017737.1_895","JAMJTW010000008.1_195"]],"hmmregion":0,
              "tree": fr"{HMM_TREE_DIR}\FlgR_db_FAMSA_gt0.1_m8ordered.tree",
              "fasta": fr"{MSA_DIR}\FlgR_db_FAMSA_gt0.1_treeordered.fasta"},
    "FliF2" : {"boundary":[["NZ_AP021861.1_3272","JABAAH010000011.1_11"]],"hmmregion":0,
              "tree": fr"{HMM_TREE_DIR}\FliF2_db_FAMSA_gt0.1_m8ordered.tree",
              "fasta": fr"{MSA_DIR}\FliF2_db_FAMSA_gt0.1_treeordered.fasta"},
    "MotB2" : {"boundary":[["VGZQ01000054.1_9","NZ_JPVZ01000003.1_410"],["VFQY01000009.1_14","DRJZ01000039.1_8",-1]],"hmmregion":0,
              "tree": fr"{HMM_TREE_DIR}\MotB2_db_FAMSA_gt0.1_m8ordered.tree",
              "fasta": fr"{MSA_DIR}\MotB2_db_FAMSA_gt0.1_treeordered.fasta"},
    "FliH2" : {"boundary":[["NZ_CP022998.1_1258","JAFLCG010000011.1_137"]],"hmmregion":0,
              "tree": fr"{HMM_TREE_DIR}\FliH2_db_FAMSA_gt0.1_m8ordered.tree",
              "fasta": fr"{MSA_DIR}\FliH2_db_FAMSA_gt0.1_treeordered.fasta"}
    
}
print(gene_boundaries.keys())

#%% Extract orthologs for every gene in `gene_boundaries`
for gene_name in gene_boundaries:
    print("Processing:", gene_name)
    gene_data = gene_boundaries[gene_name]
    hmmregion = gene_data["hmmregion"]
    boundary = gene_data["boundary"]

    if hmmregion == 1:
        hmmregion_txt = "hmmregions_"
    elif hmmregion == 0:
        hmmregion_txt = ""
    else:
        continue

    if "tree" in gene_data:
        fasta_file = gene_data["fasta"]
        tree = gene_data["tree"]
    else:
        tree = os.path.join(
            HMM_TREE_DIR,
            f"{gene_name}_hmm_E1000_db_{hmmregion_txt}FAMSA_gt0.1_hmmordered.tree",
        )
        fasta_file = os.path.join(
            MSA_DIR,
            f"{gene_name}_hmm_E1000_db_{hmmregion_txt}FAMSA_gt0.1_treeordered.fasta",
        )

    out_file = os.path.join(
        ORTHOLOG_OUT_DIR,
        f"{gene_name}_hmm_E1000_db_FAMSA_gt0.1_treeordered_orthologs.fasta",
    )
    ortholog_ids = find_leaf_nodes_within_boundaries(boundary, tree)
    filter_fasta_by_ids(fasta_file, ortholog_ids, out_file)

    is_m8_tree = "m8" in os.path.basename(tree).lower()
    if is_m8_tree:
        out_tree = derive_m8_ortholog_tree_path(gene_name, tree, HMM_TREE_DIR)
        prune_tree_by_ids(tree, ortholog_ids, out_tree)
    else:
        for hmm_tree in derive_hmm_tree_paths(gene_name, HMM_TREE_DIR):
            if not os.path.exists(hmm_tree):
                print(f"HMM tree not found for {gene_name}: {hmm_tree}")
                continue
            out_tree = hmm_tree.replace(".tree", "_orthologs.tree")
            prune_tree_by_ids(hmm_tree, ortholog_ids, out_tree)

