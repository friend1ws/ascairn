#! /usr/bin/env python

import math
from scipy.stats import nbinom
import polars as pl


def calc_loglikelihood(PR1, PR2, D_count, marker_list):

    PR12 = PR1.filter(pl.col("Marker").is_in(marker_list)) \
        .join(PR2, on="Marker", suffix = "_PR2") \
        .with_columns([
            (pl.col("prob_0") * pl.col("prob_0_PR2")).alias("mprob_0"),
            (pl.col("prob_1") * pl.col("prob_0_PR2") + pl.col("prob_0") * pl.col("prob_1_PR2")).alias("mprob_1"),
            (pl.col("prob_2") * pl.col("prob_0_PR2") + pl.col("prob_0") * pl.col("prob_2_PR2") + pl.col("prob_1") * pl.col("prob_1_PR2")).alias("mprob_2"),
            (pl.col("prob_2") * pl.col("prob_1_PR2") + pl.col("prob_1") * pl.col("prob_2_PR2")).alias("mprob_3"),
            (pl.col("prob_2") * pl.col("prob_2_PR2")).alias("mprob_4")
        ])

    PR12_count = (
        PR12.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_0") * pl.col("dprob_0") +
             pl.col("mprob_1") * pl.col("dprob_1") +
             pl.col("mprob_2") * pl.col("dprob_2") +
             pl.col("mprob_3") * pl.col("dprob_3") +
             pl.col("mprob_4") * pl.col("dprob_4")).alias("tprob_01234")
        ])
    )

    tL = PR12_count["tprob_01234"].log().sum()

    return(tL)


def calc_loglikelihood_single(PR1, D_count, marker_list):

    PR1 = PR1.filter(pl.col("Marker").is_in(marker_list)) \
        .with_columns([
            (pl.col("prob_0")).alias("mprob_0"),
            (pl.col("prob_1")).alias("mprob_1"),
            (pl.col("prob_2")).alias("mprob_2")
        ])

    PR1_count = (
        PR1.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_0") * pl.col("dprob_0") +
             pl.col("mprob_1") * pl.col("dprob_1") +
             pl.col("mprob_2") * pl.col("dprob_2")).alias("tprob_012")
        ])
    )

    tL = PR1_count["tprob_012"].log().sum()

    return(tL)


def calc_posterior_prob(PR1, PR2, D_count):

    PR12 = PR1.join(PR2, on="Marker", suffix = "_PR2") \
        .with_columns([
            (pl.col("prob_0") * pl.col("prob_0_PR2")).alias("mprob_00"),
            (pl.col("prob_0") * pl.col("prob_1_PR2")).alias("mprob_01"),
            (pl.col("prob_1") * pl.col("prob_0_PR2")).alias("mprob_10"),
            (pl.col("prob_1") * pl.col("prob_1_PR2")).alias("mprob_11"),
            (pl.col("prob_0") * pl.col("prob_2_PR2")).alias("mprob_02"),
            (pl.col("prob_2") * pl.col("prob_0_PR2")).alias("mprob_20"),
            (pl.col("prob_1") * pl.col("prob_2_PR2")).alias("mprob_12"),
            (pl.col("prob_2") * pl.col("prob_1_PR2")).alias("mprob_21"),
            (pl.col("prob_2") * pl.col("prob_2_PR2")).alias("mprob_22")
        ])

    PR12_count = PR12.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_00") * pl.col("dprob_0")).alias("tprob_00"),
            (pl.col("mprob_01") * pl.col("dprob_1")).alias("tprob_01"), 
            (pl.col("mprob_10") * pl.col("dprob_1")).alias("tprob_10"), 
            (pl.col("mprob_11") * pl.col("dprob_2")).alias("tprob_11"), 
            (pl.col("mprob_20") * pl.col("dprob_2")).alias("tprob_20"), 
            (pl.col("mprob_02") * pl.col("dprob_2")).alias("tprob_02"), 
            (pl.col("mprob_12") * pl.col("dprob_3")).alias("tprob_12"),
            (pl.col("mprob_21") * pl.col("dprob_3")).alias("tprob_21"),
            (pl.col("mprob_22") * pl.col("dprob_4")).alias("tprob_22")
        ]) \
        .with_columns([
            (pl.col("tprob_00") + pl.col("tprob_01") + pl.col("tprob_10") + pl.col("tprob_11") + \
                pl.col("tprob_20") + pl.col("tprob_02") + pl.col("tprob_12") + pl.col("tprob_21") + pl.col("tprob_22")).alias("tprob_sum")]) \
        .with_columns([
            (pl.col("tprob_00") / pl.col("tprob_sum")).alias("Prob_00"),
            (pl.col("tprob_01") / pl.col("tprob_sum")).alias("Prob_01"),
            (pl.col("tprob_10") / pl.col("tprob_sum")).alias("Prob_10"),
            (pl.col("tprob_11") / pl.col("tprob_sum")).alias("Prob_11"),
            (pl.col("tprob_02") / pl.col("tprob_sum")).alias("Prob_02"),
            (pl.col("tprob_20") / pl.col("tprob_sum")).alias("Prob_20"),
            (pl.col("tprob_12") / pl.col("tprob_sum")).alias("Prob_12"),
            (pl.col("tprob_21") / pl.col("tprob_sum")).alias("Prob_21"),
            (pl.col("tprob_22") / pl.col("tprob_sum")).alias("Prob_22")
        ]) 

    return(PR12_count)


def calc_posterior_prob_single(PR1, D_count):

    PR1 = PR1 \
        .with_columns([
            (pl.col("prob_0")).alias("mprob_0"),
            (pl.col("prob_1")).alias("mprob_1"),
            (pl.col("prob_2")).alias("mprob_2")
        ])

    PR1_count = PR1.join(D_count, on="Marker") \
        .with_columns([
            (pl.col("mprob_0") * pl.col("dprob_0")).alias("tprob_0"),
            (pl.col("mprob_1") * pl.col("dprob_1")).alias("tprob_1"), 
            (pl.col("mprob_2") * pl.col("dprob_2")).alias("tprob_2")
        ]) \
        .with_columns([
            (pl.col("tprob_0") + pl.col("tprob_1") + pl.col("tprob_2")).alias("tprob_sum")]) \
        .with_columns([
            (pl.col("tprob_0") / pl.col("tprob_sum")).alias("Prob_0"),
            (pl.col("tprob_1") / pl.col("tprob_sum")).alias("Prob_1"),
            (pl.col("tprob_2") / pl.col("tprob_sum")).alias("Prob_2")
        ]) 

    return(PR1_count)

def match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
    cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4):

    max_depth_thres = math.ceil(depth * 3)
    # max_depth_thres = 100
   
    prob_0 = [nbinom.pmf(x, nbinom_size_0, nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)) for x in range(max_depth_thres + 1)]
    prob_1 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 1.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]
    prob_2 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 2.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]
    prob_3 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 3.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]
    prob_4 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 4.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]


    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count = pl.DataFrame({
        "Marker": count["Marker"],
        "dprob_0": [prob_0[c] for c in count["Count"]],
        "dprob_1": [prob_1[c] for c in count["Count"]],
        "dprob_2": [prob_2[c] for c in count["Count"]],
        "dprob_3": [prob_3[c] for c in count["Count"]],
        "dprob_4": [prob_4[c] for c in count["Count"]]
    })

    cluster_marker_count_df = pl.read_csv(cluster_kmer_count_file, separator = '\t')


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .with_columns([pl.col("Rel_pos_std").cast(pl.Float64)]) \
            .select("Marker").unique()
            # .filter((pl.col("Rel_pos_mean") > -0.1) & (pl.col("Rel_pos_mean") < 1.1) & (pl.col("Rel_pos_std") < 0.5)) \
            # .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns([
        ((pl.col("Count_0") + pseudo_count) / 
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_0"),
        ((pl.col("Count_1") + pseudo_count) / 
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_1"),
        ((pl.col("Count_2") + pseudo_count) / 
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_2")
    ])

    # cluster_marker_count_df_2.write_csv("cluster_marker_count_df_2.tsv", separator = '\t')

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
            # print(f"({i}, {j}, {tL})")

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Cluster2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    # for debug
    # PR1_max.write_csv("PR1_max.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max.tsv", separator = '\t')

    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count) \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num", "Rel_pos_mean", "Rel_pos_std",
                    "Prob_00", "Prob_01", "Prob_10", "Prob_11", "Prob_02", "Prob_20", "Prob_12", "Prob_21", "Prob_22"])


    PR12_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    if cluster_haplotype_file is None: return


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

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns([
        # Creating binary columns based on `Count` values
            (pl.col("Count") == 0).cast(pl.Int64).alias("Count_0"),
            (pl.col("Count") == 1).cast(pl.Int64).alias("Count_1"),
            (pl.col("Count") == 2).cast(pl.Int64).alias("Count_2"),
        ])
        .drop("Count")
        .with_columns([
            (pl.col("Count_0") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_0"),
            (pl.col("Count_1") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_1"),
            (pl.col("Count_2") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_2"),
        ])

    )


    target_cluster_1 = D_LL[0, "Cluster1"]
    target_cluster_2 = D_LL[0, "Cluster2"]

    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select([
            "Marker",
            pl.col("prob_0").alias("prob_0_c"),
            pl.col("prob_1").alias("prob_1_c"),
            pl.col("prob_2").alias("prob_2_c")
        ])

    PR2_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select([
            "Marker",
            pl.col("prob_0").alias("prob_0_c"),
            pl.col("prob_1").alias("prob_1_c"),
            pl.col("prob_2").alias("prob_2_c")
        ])  

    # marker_list.to_frame().write_csv("marker_list.tsv", separator = '\t')
    # PR2_c.write_csv("PR1_c.tsv", separator = '\t')
    # hap_marker_count_df_2.write_csv("hap_marker_count_df_2.tsv", separator = '\t')

    cluster_haplotype_list = pl.read_csv(cluster_haplotype_file, separator = '\t')

    target_cluster_hap_1 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    target_cluster_hap_2 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_2) \
        .select("Haplotype") \
        .to_series() \
        .to_list()


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
                    .select([
                        "Marker",
                        pl.col("prob_0").alias("prob_0_h"),
                        pl.col("prob_1").alias("prob_1_h"),
                        pl.col("prob_2").alias("prob_2_h")
                    ]) \
                    .join(PR1_c, on="Marker", how="inner") \
                    .with_columns([
                        ((1 - cluster_ratio) * pl.col("prob_0_h") + cluster_ratio * pl.col("prob_0_c")).alias("prob_0"),
                        ((1 - cluster_ratio) * pl.col("prob_1_h") + cluster_ratio * pl.col("prob_1_c")).alias("prob_1"),
                        ((1 - cluster_ratio) * pl.col("prob_2_h") + cluster_ratio * pl.col("prob_2_c")).alias("prob_2")
                    ])

            else:
                PR1 = PR1_c \
                    .select([
                        "Marker",
                        pl.col("prob_0_c").alias("prob_0"),
                        pl.col("prob_1_c").alias("prob_1"),
                        pl.col("prob_2_c").alias("prob_2")
                    ])

            if target_cluster_hap_2[j] not in [None, "NA", "None", "NONE"]:

                PR2 = hap_marker_count_df_2 \
                    .filter(pl.col("Haplotype") == target_cluster_hap_2[j]) \
                    .select([
                        "Marker",
                        pl.col("prob_0").alias("prob_0_h"),
                        pl.col("prob_1").alias("prob_1_h"),
                        pl.col("prob_2").alias("prob_2_h")
                    ]) \
                    .join(PR2_c, on="Marker", how="inner") \
                    .with_columns([
                        ((1 - cluster_ratio) * pl.col("prob_0_h") + cluster_ratio * pl.col("prob_0_c")).alias("prob_0"),
                        ((1 - cluster_ratio) * pl.col("prob_1_h") + cluster_ratio * pl.col("prob_1_c")).alias("prob_1"),
                        ((1 - cluster_ratio) * pl.col("prob_2_h") + cluster_ratio * pl.col("prob_2_c")).alias("prob_2")
                    ])

            else:
                PR2 = PR2_c \
                    .select([
                        "Marker",
                        pl.col("prob_0_c").alias("prob_0"),
                        pl.col("prob_1_c").alias("prob_1"),
                        pl.col("prob_2_c").alias("prob_2")
                    ])

     
            tL = calc_loglikelihood(PR1, PR2, D_count, marker_list)
            cl1_list.append(target_cluster_hap_1[i])
            cl2_list.append(target_cluster_hap_2[j])
            LL_list.append(tL)
            # print(f"({target_cluster_hap_1[i]}, {target_cluster_hap_2[j]}, {tL})")

            if tL > LL_max:
                PR1_max = PR1
                PR2_max = PR2
                LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Haplotype2": cl2_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')


    hap_info = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"], D_LL2[0, "Haplotype2"]])) \
        .group_by(["Marker", "Haplotype", "Marker_num"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"]) 
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    hap_info2 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype2"]) 
        .rename({"Haplotype": "Haplotype2", "Mean_marker_pos": "Mean_marker_pos2", "Marker_count": "Marker_count2"})
    )



    PR12_count = calc_posterior_prob(PR1_max, PR2_max, D_count) \
        .join(hap_info1, on = "Marker", how = "left") \
        .join(hap_info2, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA"), 
                 pl.col("Haplotype2").fill_null("NA"), pl.col("Mean_marker_pos2").fill_null("NA"), pl.col("Marker_count2").fill_null("NA"),
                 "Prob_00", "Prob_01", "Prob_10", "Prob_11", "Prob_02", "Prob_20", "Prob_12", "Prob_21", "Prob_22"])

    # PR1_max.join(PR2_max, on="Marker", suffix = "_PR2").write_csv("PR12_max.tsv", separator = '\t')
    # PR1_max.write_csv("PR1_max_hap.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max_hap.tsv", separator = '\t')
    PR12_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')



def match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
    cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4):

    max_depth_thres = math.ceil(depth * 1.5 * 0.5) 
    # max_depth_thres = 200

    prob_0 = [nbinom.pmf(x, nbinom_size_0, nbinom_size_0 / (nbinom_size_0 + nbinom_mu_0)) for x in range(max_depth_thres + 1)]
    prob_1 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 1.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]
    prob_2 = [nbinom.pmf(x, nbinom_size, nbinom_size / (nbinom_size + 2.0 * nbinom_mu_unit * depth)) for x in range(max_depth_thres + 1)]


    count = pl.read_csv(kmer_count_file, separator = '\t', new_columns = ["Marker", "Count"]) \
            .filter(pl.col("Count") <= max_depth_thres)

    D_count = pl.DataFrame({
        "Marker": count["Marker"],
        "dprob_0": [prob_0[c] for c in count["Count"]],
        "dprob_1": [prob_1[c] for c in count["Count"]],
        "dprob_2": [prob_2[c] for c in count["Count"]],
    })

    cluster_marker_count_df = pl.read_csv(cluster_kmer_count_file, separator = '\t')


    marker_list =  pl.Series(cluster_marker_count_df \
            .filter((pl.col("Marker_num") >= 2) & (pl.col("Hap_minus_marker_num") >= 2)) \
            .with_columns([pl.col("Rel_pos_std").cast(pl.Float64)]) \
            .filter((pl.col("Rel_pos_mean") > -0.1) & (pl.col("Rel_pos_mean") < 1.1) & (pl.col("Rel_pos_std") < 0.5)) \
            .select("Marker").unique()
        )

    cluster_num = cluster_marker_count_df["Cluster"].max()

    cluster_marker_count_df_2 = cluster_marker_count_df.with_columns([
        ((pl.col("Count_0") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_0"),
        ((pl.col("Count_1") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_1"),
        ((pl.col("Count_2") + pseudo_count) /
         (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2") + 3 * pseudo_count)).alias("prob_2")
    ])

    # cluster_marker_count_df_2.write_csv("cluster_marker_count_df_2.tsv", separator = '\t')

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(1, cluster_num + 1):

        PR1 = cluster_marker_count_df_2.filter(pl.col("Cluster") == i)

        tL = calc_loglikelihood_single(PR1, D_count, marker_list)
        cl1_list.append(i)
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL = pl.DataFrame({"Cluster1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)
    D_LL.write_csv(output_prefix + ".cluster.hap_pair.txt", separator = '\t')

    PR1_count = calc_posterior_prob_single(PR1_max, D_count) \
                    .select(["Marker", "Marker_num", "Hap_minus_marker_num", "Rel_pos_mean", "Rel_pos_std",
                    "Prob_0", "Prob_1", "Prob_2"])

    PR1_count.write_csv(output_prefix + ".cluster.marker_prob.txt", separator = '\t')


    if cluster_haplotype_file is None: return


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

    hap_marker_count_df_2 = (
        all_combinations.join(hap_marker_count_df_2_tmp, on=["Marker", "Haplotype"], how="left")
        .select(["Marker", "Haplotype", pl.col("Count").fill_null(0)])
        .with_columns([
        # Creating binary columns based on `Count` values
            (pl.col("Count") == 0).cast(pl.Int64).alias("Count_0"),
            (pl.col("Count") == 1).cast(pl.Int64).alias("Count_1"),
            (pl.col("Count") == 2).cast(pl.Int64).alias("Count_2"),
        ])
        .drop("Count")
        .with_columns([
            (pl.col("Count_0") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_0"),
            (pl.col("Count_1") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_1"),
            (pl.col("Count_2") / (pl.col("Count_0") + pl.col("Count_1") + pl.col("Count_2"))).alias("prob_2"),
        ])

    )


    target_cluster_1 = D_LL[0, "Cluster1"]

    PR1_c = cluster_marker_count_df_2 \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select([
            "Marker",
            pl.col("prob_0").alias("prob_0_c"),
            pl.col("prob_1").alias("prob_1_c"),
            pl.col("prob_2").alias("prob_2_c")
        ])

    cluster_haplotype_list = pl.read_csv(cluster_haplotype_file, separator = '\t')

    target_cluster_hap_1 = cluster_haplotype_list \
        .filter(pl.col("Cluster") == target_cluster_1) \
        .select("Haplotype") \
        .to_series() \
        .to_list()

    cl1_list = []
    LL_list = []

    LL_max = -float("Inf")
    PR1_max = None


    for i in range(len(target_cluster_hap_1)):

        if target_cluster_hap_1[i] not in [None, "NA", "None", "NONE"]:

            PR1 = hap_marker_count_df_2 \
                .filter(pl.col("Haplotype") == target_cluster_hap_1[i]) \
                .select([
                    "Marker",
                    pl.col("prob_0").alias("prob_0_h"),
                    pl.col("prob_1").alias("prob_1_h"),
                    pl.col("prob_2").alias("prob_2_h")
                ]) \
                .join(PR1_c, on="Marker", how="inner") \
                .with_columns([
                    ((1 - cluster_ratio) * pl.col("prob_0_h") + cluster_ratio * pl.col("prob_0_c")).alias("prob_0"),
                    ((1 - cluster_ratio) * pl.col("prob_1_h") + cluster_ratio * pl.col("prob_1_c")).alias("prob_1"),
                    ((1 - cluster_ratio) * pl.col("prob_2_h") + cluster_ratio * pl.col("prob_2_c")).alias("prob_2")
                ])

        else:
            PR1 = PR1_c \
                .select([
                    "Marker",
                    pl.col("prob_0_c").alias("prob_0"),
                    pl.col("prob_1_c").alias("prob_1"),
                    pl.col("prob_2_c").alias("prob_2")
                ])

        tL = calc_loglikelihood_single(PR1, D_count, marker_list)
        cl1_list.append(target_cluster_hap_1[i])
        LL_list.append(tL)

        if tL > LL_max:
            PR1_max = PR1
            LL_max = tL


    D_LL2 = pl.DataFrame({"Haplotype1": cl1_list, "Loglikelihood": LL_list}).sort("Loglikelihood", descending = True)

    D_LL2.write_csv(output_prefix + ".haplotype.hap_pair.txt", separator = '\t')



    hap_info = pl.read_csv(kmer_info_file, separator = '\t', infer_schema_length = None) \
        .filter(pl.col("Haplotype").is_in([D_LL2[0, "Haplotype1"]])) \
        .group_by(["Marker", "Haplotype", "Marker_num"]) \
        .agg([
            pl.col("Marker_pos").mean().alias("Mean_marker_pos"),
            pl.col("Marker_pos").count().alias("Marker_count")
        ])

    hap_info1 = (
        hap_info.filter(pl.col("Haplotype") == D_LL2[0, "Haplotype1"])
        .rename({"Haplotype": "Haplotype1", "Mean_marker_pos": "Mean_marker_pos1", "Marker_count": "Marker_count1"})
    )

    PR1_count = calc_posterior_prob_single(PR1_max, D_count) \
        .join(hap_info1, on = "Marker", how = "left") \
        .select(["Marker", pl.col("Haplotype1").fill_null("NA"), pl.col("Mean_marker_pos1").fill_null("NA"), pl.col("Marker_count1").fill_null("NA"),
                 "Prob_0", "Prob_1", "Prob_2"])

    # PR1_max.join(PR2_max, on="Marker", suffix = "_PR2").write_csv("PR12_max.tsv", separator = '\t')
    # PR1_max.write_csv("PR1_max_hap.tsv", separator = '\t')
    # PR2_max.write_csv("PR2_max_hap.tsv", separator = '\t')
    PR1_count.write_csv(output_prefix + ".haplotype.marker_prob.txt", separator = '\t')


if __name__ == "__main__":

    # print("test")

    kmer_count_file = "../out/out.kmer_count.txt"
    output_prefix = "output_pl/out.match/chr1"
    kmer_info_file = "../ascairn/data/kmer_info/chr1.kmer_info.txt.gz"
    cluster_kmer_count_file = "../ascairn/data/cluster/chr1.cluster_marker_count.txt.gz"
    depth = 41.68
    cluster_haplotype_file = "../ascairn/data/cluster/chr1.hap_cluster.txt"

    match_cluster_haplotype(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
        cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)


    kmer_count_file = "../out/out.kmer_count.txt"
    output_prefix = "output_pl/out.match/chrX"
    kmer_info_file = "../ascairn/data/kmer_info/chrX.kmer_info.txt.gz"
    cluster_kmer_count_file = "../ascairn/data/cluster/chrX.cluster_marker_count.txt.gz"
    depth = 41.68 
    cluster_haplotype_file = "../ascairn/data/cluster/chrX.hap_cluster.txt"
 
    match_cluster_haplotype_single(kmer_count_file, output_prefix, kmer_info_file, cluster_kmer_count_file, depth,
        cluster_haplotype_file, cluster_ratio = 0.1, pseudo_count = 0.1, nbinom_size_0 = 0.5, nbinom_size = 8, nbinom_mu_0 = 0.8, nbinom_mu_unit = 0.4)



