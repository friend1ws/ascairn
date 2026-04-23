#! /usr/bin/env python

import math
from scipy.stats import nbinom
import polars as pl


def calc_loglikelihood(PR1, PR2, D_count, marker_list, max_copy_number=2):
    # mprob_k = sum of prob_i * prob_j where i + j = k
    mprob_exprs = []
    for k in range(2 * max_copy_number + 1):
        terms = [pl.col(f"prob_{i}") * pl.col(f"prob_{k-i}_PR2")
                 for i in range(max_copy_number + 1) if 0 <= k - i <= max_copy_number]
        mprob_exprs.append(sum(terms).alias(f"mprob_{k}"))

    PR12 = PR1.filter(pl.col("Marker").is_in(marker_list)) \
        .join(PR2, on="Marker", suffix="_PR2") \
        .with_columns(mprob_exprs)

    tprob_expr = sum([pl.col(f"mprob_{k}") * pl.col(f"dprob_{k}")
                      for k in range(2 * max_copy_number + 1)]).alias("tprob")

    PR12_count = PR12.join(D_count, on="Marker").with_columns([tprob_expr])

    tL = PR12_count["tprob"].log().sum()

    return(tL)


def calc_loglikelihood_single(PR1, D_count, marker_list, max_copy_number=2):

    PR1 = PR1.filter(pl.col("Marker").is_in(marker_list))

    tprob_expr = sum([pl.col(f"prob_{k}") * pl.col(f"dprob_{k}")
                      for k in range(max_copy_number + 1)]).alias("tprob")

    PR1_count = PR1.join(D_count, on="Marker").with_columns([tprob_expr])

    tL = PR1_count["tprob"].log().sum()

    return(tL)


def calc_posterior_prob(PR1, PR2, D_count, max_copy_number=2):
    # Generate (N+1)² mprob states: mprob_ij = prob_i * prob_j
    mprob_exprs = [(pl.col(f"prob_{i}") * pl.col(f"prob_{j}_PR2")).alias(f"mprob_{i}{j}")
                   for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    PR12 = PR1.join(PR2, on="Marker", suffix="_PR2").with_columns(mprob_exprs)

    # tprob_ij = mprob_ij * dprob_{i+j}
    tprob_exprs = [(pl.col(f"mprob_{i}{j}") * pl.col(f"dprob_{i+j}")).alias(f"tprob_{i}{j}")
                   for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    # Sum all tprob values
    tprob_sum_expr = sum([pl.col(f"tprob_{i}{j}") for i in range(max_copy_number + 1)
                          for j in range(max_copy_number + 1)]).alias("tprob_sum")

    # Normalize to get posterior probabilities
    prob_exprs = [(pl.col(f"tprob_{i}{j}") / pl.col("tprob_sum")).alias(f"Prob_{i}{j}")
                  for i in range(max_copy_number + 1) for j in range(max_copy_number + 1)]

    PR12_count = PR12.join(D_count, on="Marker") \
        .with_columns(tprob_exprs) \
        .with_columns([tprob_sum_expr]) \
        .with_columns(prob_exprs)

    return(PR12_count)


def calc_posterior_prob_single(PR1, D_count, max_copy_number=2):

    tprob_exprs = [(pl.col(f"prob_{k}") * pl.col(f"dprob_{k}")).alias(f"tprob_{k}")
                   for k in range(max_copy_number + 1)]

    tprob_sum_expr = sum([pl.col(f"tprob_{k}") for k in range(max_copy_number + 1)]).alias("tprob_sum")

    prob_exprs = [(pl.col(f"tprob_{k}") / pl.col("tprob_sum")).alias(f"Prob_{k}")
                  for k in range(max_copy_number + 1)]

    PR1_count = PR1.join(D_count, on="Marker") \
        .with_columns(tprob_exprs) \
        .with_columns([tprob_sum_expr]) \
        .with_columns(prob_exprs)

    return(PR1_count)


def calc_rel_pos(kmer_info_file):
    """Calculate Rel_pos_mean and Rel_pos_std per Marker from kmer_info."""
    kmer_info_df = pl.read_csv(kmer_info_file, infer_schema_length=None, separator='\t')
    rel_pos_df = kmer_info_df.with_columns(
        (pl.col("Marker_pos") / pl.col("Contig_len")).alias("Rel_pos")
    ).group_by("Marker").agg([
        pl.col("Rel_pos").mean().alias("Rel_pos_mean"),
        pl.col("Rel_pos").std().alias("Rel_pos_std"),
    ])
    return rel_pos_df


def build_cluster_marker_count(kmer_info_file, hap_info_file):
    """Build cluster_marker_count DataFrame on-the-fly from kmer_info and hap_info."""

    kmer_info_df = pl.read_csv(kmer_info_file, infer_schema_length=None, separator='\t')
    hap_info_df = pl.read_csv(hap_info_file, separator='\t')

    # Count occurrences per Marker x Haplotype (copy number)
    hap_marker_count = kmer_info_df.select("Marker", "Haplotype") \
        .group_by("Marker", "Haplotype") \
        .agg(pl.len().alias("copy_number"))

    # Auto-detect max_copy_number from actual data
    max_copy_number = hap_marker_count["copy_number"].max()
    if max_copy_number >= 3:
        raise ValueError(f"Max copy number in kmer_info is {max_copy_number} (>= 3). This is not supported.")

    # Get all unique markers and haplotypes for complete combinations
    unique_markers = hap_marker_count.select("Marker").unique()
    unique_haplotypes = hap_info_df.select("Haplotype").unique()
    all_combinations = unique_markers.join(unique_haplotypes, how="cross")

    # Fill missing with copy_number=0
    hap_marker_full = all_combinations.join(hap_marker_count, on=["Marker", "Haplotype"], how="left") \
        .with_columns(pl.col("copy_number").fill_null(0))

    # Join with cluster info
    hap_marker_full = hap_marker_full.join(hap_info_df.select("Haplotype", "Cluster"), on="Haplotype")

    # Count per Cluster x Marker x copy_number
    cluster_marker_count = hap_marker_full.group_by("Cluster", "Marker", "copy_number") \
        .agg(pl.len().alias("Count"))

    # Pivot to Count_0, Count_1, ..., Count_N
    cluster_marker_count = cluster_marker_count.pivot(
        on="copy_number", index=["Cluster", "Marker"], values="Count"
    )

    # Rename pivoted columns and fill nulls
    for i in range(max_copy_number + 1):
        col_name = str(i)
        target_name = f"Count_{i}"
        if col_name in cluster_marker_count.columns:
            cluster_marker_count = cluster_marker_count.rename({col_name: target_name})
        else:
            cluster_marker_count = cluster_marker_count.with_columns(pl.lit(0).alias(target_name))
        cluster_marker_count = cluster_marker_count.with_columns(pl.col(target_name).fill_null(0))

    # Add Marker_num and Hap_minus_marker_num
    # Marker_num: number of distinct haplotypes that have the marker (copy_number >= 1)
    # Computed from the unique Marker x Haplotype pairs in hap_marker_count (already deduplicated)
    total_hap = hap_info_df.shape[0]
    marker_num_df = hap_marker_count.group_by("Marker") \
        .agg(pl.len().alias("Marker_num"))
    marker_num_df = marker_num_df.with_columns(
        (pl.lit(total_hap) - pl.col("Marker_num")).alias("Hap_minus_marker_num")
    )

    cluster_marker_count = cluster_marker_count.join(marker_num_df, on="Marker", how="left") \
        .with_columns([
            pl.col("Marker_num").fill_null(0),
            pl.col("Hap_minus_marker_num").fill_null(total_hap)
        ])

    return cluster_marker_count, max_copy_number


def match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, hap_info_file, depth,
    cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0_unit = 0.8 / 30, nbinom_mu_unit = 0.4,
    hap_candidates_file = None):

    cluster_marker_count_df, max_copy_number = build_cluster_marker_count(kmer_info_file, hap_info_file)

    max_depth_thres = math.ceil(depth * 3)

    # prob_0: background distribution
    nbinom_mu_0 = nbinom_mu_0_unit * depth
    nbinom_p_0 = nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)
    prob_list = [[nbinom.pmf(x, nbinom_size_0, nbinom_p_0) for x in range(max_depth_thres + 1)]]
    # prob_1 through prob_{2*max_copy_number}: signal distributions
    for i in range(1, 2 * max_copy_number + 1):
        nbinom_p = nbinom_size / (nbinom_size + i * nbinom_mu_unit * depth)
        prob_list.append([nbinom.pmf(x, nbinom_size, nbinom_p) for x in range(max_depth_thres + 1)])

    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count_dict = {"Marker": count["Marker"]}
    for i, prob in enumerate(prob_list):
        D_count_dict[f"dprob_{i}"] = [prob[c] for c in count["Count"]]
    D_count = pl.DataFrame(D_count_dict)


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    total = sum(count_cols) + (max_copy_number + 1) * pseudo_count
    prob_exprs = [((pl.col(f"Count_{i}") + pseudo_count) / total).alias(f"prob_{i}")
                  for i in range(max_copy_number + 1)]
    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns(prob_exprs)

    cl1_list = []
    cl2_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None
    PR2_max = None

    for i in range(1, cluster_num + 1):
        for j in range(i, cluster_num + 1):

            PR1 = cluster_marker_count_df_2.filter(pl.col("Cluster") == i)
            PR2 = cluster_marker_count_df_2.filter(pl.col("Cluster") == j)

            tL = calc_loglikelihood(PR1, PR2, D_count, marker_list)
            cl1_list.append(i)
            cl2_list.append(j)
            LL_list.append(tL)

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Cluster2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    rel_pos_df = calc_rel_pos(kmer_info_file)

    prob_cols = [f"Prob_{i}{j}" for i in range(max_copy_number + 1)
                 for j in range(max_copy_number + 1)]
    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count, max_copy_number) \
                    .join(rel_pos_df, on="Marker", how="left") \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num", "Rel_pos_mean", "Rel_pos_std"] + prob_cols)


    PR12_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    hap_marker_count_df_2_tmp = pl.read_csv(kmer_info_file, infer_schema_length = None, separator = '\t') \
        .select("Marker", "Haplotype") \
        .group_by("Marker", "Haplotype") \
        .agg(pl.len().alias("Count"))

    unique_markers = hap_marker_count_df_2_tmp.select("Marker").unique()
    unique_haplotypes = hap_marker_count_df_2_tmp.select("Haplotype").unique()

    # Create all combinations of Markers and Haplotypes
    all_combinations = (
        unique_markers.join(unique_haplotypes, how="cross")
    )

    # Count_0 through Count_N: exact match
    count_exprs = [(pl.col("Count") == i).cast(pl.Int64).alias(f"Count_{i}")
                   for i in range(max_copy_number + 1)]

    hap_count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    hap_total = sum(hap_count_cols)
    hap_prob_exprs = [(pl.col(f"Count_{i}") / hap_total).alias(f"prob_{i}")
                      for i in range(max_copy_number + 1)]

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns(count_exprs)
        .drop("Count")
        .with_columns(hap_prob_exprs)
    )


    target_cluster_1 = D_LL[0, "Cluster1"]
    target_cluster_2 = D_LL[0, "Cluster2"]

    prob_c_cols = ["Marker"] + [pl.col(f"prob_{i}").alias(f"prob_{i}_c") for i in range(max_copy_number + 1)]
    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select(prob_c_cols)

    PR2_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select(prob_c_cols)

    hap_info_df = pl.read_csv(hap_info_file, separator = '\t')

    # Load haplotype candidate list if provided
    hap_candidates = None
    if hap_candidates_file is not None:
        with open(hap_candidates_file) as f:
            hap_candidates = set(line.strip() for line in f if line.strip())

    target_cluster_hap_1 = hap_info_df \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    target_cluster_hap_2 = hap_info_df \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    if hap_candidates is not None:
        target_cluster_hap_1 = [h for h in target_cluster_hap_1 if h in hap_candidates]
        target_cluster_hap_2 = [h for h in target_cluster_hap_2 if h in hap_candidates]

    cl1_list = []
    cl2_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None
    PR2_max = None


    for i in range(len(target_cluster_hap_1)):
        for j in range(len(target_cluster_hap_2)):

            if target_cluster_hap_1[i] not in [None, "NA", "None", "NONE"]:

                PR1 = hap_marker_count_df_2 \
                    .filter(pl.col("Haplotype") == target_cluster_hap_1[i]) \
                    .select(["Marker"] + [pl.col(f"prob_{k}").alias(f"prob_{k}_h") for k in range(max_copy_number + 1)]) \
                    .join(PR1_c, on="Marker", how="inner") \
                    .with_columns([((1 - cluster_ratio) * pl.col(f"prob_{k}_h") + cluster_ratio * pl.col(f"prob_{k}_c")).alias(f"prob_{k}")
                                   for k in range(max_copy_number + 1)])

            else:
                PR1 = PR1_c \
                    .select(["Marker"] + [pl.col(f"prob_{k}_c").alias(f"prob_{k}") for k in range(max_copy_number + 1)])

            if target_cluster_hap_2[j] not in [None, "NA", "None", "NONE"]:

                PR2 = hap_marker_count_df_2 \
                    .filter(pl.col("Haplotype") == target_cluster_hap_2[j]) \
                    .select(["Marker"] + [pl.col(f"prob_{k}").alias(f"prob_{k}_h") for k in range(max_copy_number + 1)]) \
                    .join(PR2_c, on="Marker", how="inner") \
                    .with_columns([((1 - cluster_ratio) * pl.col(f"prob_{k}_h") + cluster_ratio * pl.col(f"prob_{k}_c")).alias(f"prob_{k}")
                                   for k in range(max_copy_number + 1)])

            else:
                PR2 = PR2_c \
                    .select(["Marker"] + [pl.col(f"prob_{k}_c").alias(f"prob_{k}") for k in range(max_copy_number + 1)])

     
            tL = calc_loglikelihood(PR1, PR2, D_count, marker_list)
            cl1_list.append(target_cluster_hap_1[i])
            cl2_list.append(target_cluster_hap_2[j])
            LL_list.append(tL)

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Haplotype2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')


    marker_hap_summary = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"], D_LL2[0, "Haplotype2"]])) \
        .group_by(["Marker", "Haplotype"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        marker_hap_summary.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"])
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    hap_info2 = (
        marker_hap_summary.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype2"])
        .rename({"Haplotype": "Haplotype2", "Mean_marker_pos": "Mean_marker_pos2", "Marker_count": "Marker_count2"})
    )



    prob_cols = [f"Prob_{i}{j}" for i in range(max_copy_number + 1)
                 for j in range(max_copy_number + 1)]
    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count, max_copy_number) \
        .join(hap_info1, on = "Marker", how = "left") \
        .join(hap_info2, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA"),
                 pl.col("Haplotype2").fill_null("NA"), pl.col("Mean_marker_pos2").fill_null("NA"), pl.col("Marker_count2").fill_null("NA")] + prob_cols)

    PR12_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')

    # Write cen_type.txt with hap_info annotations auto-expanded
    best_hap1 = D_LL2[0, "Haplotype1"]
    best_hap2 = D_LL2[0, "Haplotype2"]
    best_cl1 = D_LL[0, "Cluster1"]
    best_cl2 = D_LL[0, "Cluster2"]

    extra_cols = [c for c in hap_info_df.columns if c not in ("Haplotype", "Cluster")]
    row1 = hap_info_df.filter(pl.col("Haplotype") == best_hap1)
    row2 = hap_info_df.filter(pl.col("Haplotype") == best_hap2)

    result_data = {"Cluster_1": [best_cl1], "Cluster_2": [best_cl2],
                   "Haplotype_1": [best_hap1], "Haplotype_2": [best_hap2]}
    for col in extra_cols:
        result_data[f"{col}_1"] = [row1[0, col] if row1.height > 0 else "NA"]
        result_data[f"{col}_2"] = [row2[0, col] if row2.height > 0 else "NA"]

    pl.DataFrame(result_data).write_csv(output_prefix + ".cen_type.txt", separator='\t')


def match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, hap_info_file, depth,
    cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0_unit = 0.8 / 30, nbinom_mu_unit = 0.4,
    hap_candidates_file = None):

    cluster_marker_count_df, max_copy_number = build_cluster_marker_count(kmer_info_file, hap_info_file)

    max_depth_thres = math.ceil(depth * 1.5)

    # prob_0: background distribution
    nbinom_mu_0 = nbinom_mu_0_unit * depth
    nbinom_p_0 = nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)
    prob_list = [[nbinom.pmf(x, nbinom_size_0, nbinom_p_0) for x in range(max_depth_thres + 1)]]
    # prob_1 through prob_{max_copy_number}: signal distributions (single haplotype)
    for i in range(1, max_copy_number + 1):
        nbinom_p = nbinom_size / (nbinom_size + i * nbinom_mu_unit * depth)
        prob_list.append([nbinom.pmf(x, nbinom_size, nbinom_p) for x in range(max_depth_thres + 1)])

    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count_dict = {"Marker": count["Marker"]}
    for i, prob in enumerate(prob_list):
        D_count_dict[f"dprob_{i}"] = [prob[c] for c in count["Count"]]
    D_count = pl.DataFrame(D_count_dict)


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    total = sum(count_cols) + (max_copy_number + 1) * pseudo_count
    prob_exprs = [((pl.col(f"Count_{i}") + pseudo_count) / total).alias(f"prob_{i}")
                  for i in range(max_copy_number + 1)]
    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns(prob_exprs)

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(1, cluster_num + 1):

        PR1 = cluster_marker_count_df_2.filter(pl.col("Cluster") == i)

        tL = calc_loglikelihood_single(PR1, D_count, marker_list, max_copy_number)
        cl1_list.append(i)
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)
    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    rel_pos_df = calc_rel_pos(kmer_info_file)

    prob_cols_single = [f"Prob_{i}" for i in range(max_copy_number + 1)]
    PR1_count = calc_posterior_prob_single(PR1_max, D_count, max_copy_number) \
                    .join(rel_pos_df, on="Marker", how="left") \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num",
                    "Rel_pos_mean", "Rel_pos_std"] + prob_cols_single)

    PR1_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    hap_marker_count_df_2_tmp = pl.read_csv(kmer_info_file, infer_schema_length = None, separator = '\t') \
        .select("Marker", "Haplotype") \
        .group_by("Marker", "Haplotype") \
        .agg(pl.len().alias("Count"))

    unique_markers = hap_marker_count_df_2_tmp.select("Marker").unique()
    unique_haplotypes = hap_marker_count_df_2_tmp.select("Haplotype").unique()

    # Create all combinations of Markers and Haplotypes
    all_combinations = (
        unique_markers.join(unique_haplotypes, how="cross")
    )

    count_exprs = [(pl.col("Count") == i).cast(pl.Int64).alias(f"Count_{i}")
                   for i in range(max_copy_number + 1)]

    hap_count_cols = [pl.col(f"Count_{i}") for i in range(max_copy_number + 1)]
    hap_total = sum(hap_count_cols)
    hap_prob_exprs = [(pl.col(f"Count_{i}") / hap_total).alias(f"prob_{i}")
                      for i in range(max_copy_number + 1)]

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns(count_exprs)
        .drop("Count")
        .with_columns(hap_prob_exprs)
    )


    target_cluster_1 = D_LL[0, "Cluster1"]

    prob_c_cols = ["Marker"] + [pl.col(f"prob_{i}").alias(f"prob_{i}_c") for i in range(max_copy_number + 1)]
    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select(prob_c_cols)

    hap_info_df = pl.read_csv(hap_info_file, separator = '\t')

    # Load haplotype candidate list if provided
    hap_candidates = None
    if hap_candidates_file is not None:
        with open(hap_candidates_file) as f:
            hap_candidates = set(line.strip() for line in f if line.strip())

    target_cluster_hap_1 = hap_info_df \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    if hap_candidates is not None:
        target_cluster_hap_1 = [h for h in target_cluster_hap_1 if h in hap_candidates]

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(len(target_cluster_hap_1)):

        if target_cluster_hap_1[i] not in [None, "NA", "None", "NONE"]:

            PR1 = hap_marker_count_df_2 \
                .filter(pl.col("Haplotype") == target_cluster_hap_1[i]) \
                .select(["Marker"] + [pl.col(f"prob_{k}").alias(f"prob_{k}_h") for k in range(max_copy_number + 1)]) \
                .join(PR1_c, on="Marker", how="inner") \
                .with_columns([((1 - cluster_ratio) * pl.col(f"prob_{k}_h") + cluster_ratio * pl.col(f"prob_{k}_c")).alias(f"prob_{k}")
                               for k in range(max_copy_number + 1)])

        else:
            PR1 = PR1_c \
                .select(["Marker"] + [pl.col(f"prob_{k}_c").alias(f"prob_{k}") for k in range(max_copy_number + 1)])

        tL = calc_loglikelihood_single(PR1, D_count, marker_list, max_copy_number)
        cl1_list.append(target_cluster_hap_1[i])
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')



    marker_hap_summary = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"]])) \
        .group_by(["Marker", "Haplotype"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        marker_hap_summary.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"])
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    PR1_count = calc_posterior_prob_single(PR1_max, D_count, max_copy_number) \
        .join(hap_info1, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA")]
                 + prob_cols_single)

    PR1_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')

    # Write cen_type.txt with hap_info annotations auto-expanded
    best_hap1 = D_LL2[0, "Haplotype1"]
    best_cl1 = D_LL[0, "Cluster1"]

    extra_cols = [c for c in hap_info_df.columns if c not in ("Haplotype", "Cluster")]
    row1 = hap_info_df.filter(pl.col("Haplotype") == best_hap1)

    result_data = {"Cluster_1": [best_cl1], "Cluster_2": ["NA"],
                   "Haplotype_1": [best_hap1], "Haplotype_2": ["NA"]}
    for col in extra_cols:
        result_data[f"{col}_1"] = [row1[0, col] if row1.height > 0 else "NA"]
        result_data[f"{col}_2"] = ["NA"]

    pl.DataFrame(result_data).write_csv(output_prefix + ".cen_type.txt", separator='\t')
