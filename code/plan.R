
c_dir <- fs::path_expand(data_config$cause_dir)

eqtl_dir <- fs::path_expand(data_config$eqtl_dir)
all_files <- fs::path_file(fs::dir_ls(c_dir, glob = "*cause.RDS"))


all(file.exists(fs::path(c_dir, all_files)))
eqtl_files <- fs::path_file(fs::dir_ls(eqtl_dir, glob = data_config$eqtl_ext))
plan <- drake_plan(
    causedata = target(readRDS(fs::path(c_dir, all_files[x])),
                       transform = map(x = !!(seq_along(all_files)))),
    medz = target(median_Z(causedata), transform = map(causedata)),
    cause_df = target(scause_df(causedata$data), transform = map(causedata)),
    all_snps = target(unique_col(cause_df, col = "snp"),
                      transform = combine(cause_df)),
    snp_df = target(annotate_snp_loc(all_snps), hpc = FALSE),
    gene_df = target(retrieve_genes(), hpc = FALSE),
    anno_df = nearest_gene(snp_df, gene_df = gene_df),
    t_eqtl_df = target(read_eqtl_f(fs::path(eqtl_dir, eqtl_files[i])),
                       transform = map(i = !!seq_along(eqtl_files))),
    s_eqtl_df = target(subset_eqtl(anno_df, t_eqtl_df = t_eqtl_df),
                       transform = map(t_eqtl_df)),
    eqtl_df = target(bind_rows(s_eqtl_df),
                     transform = combine(s_eqtl_df)),
    genelist = unique_col(eqtl_df, anno_df, col = "ensembl_gene_stable_id"),
    godf = target(retrieve_go_ids(genelist), hpc = FALSE),
    snp_gene_df = finalize_snp_gene(anno_df, eqtl_df, godf),
    sgm = snp_gene_mat(snp_gene_df, all_snps),
    res_df = target(
        enrich_fun(cause_df = cause_df, mz = medz, snp_gene_mat = sgm, file = all_files[x]),
        transform = map(cause_df, medz, x = !!(seq_along(all_files)))
    ),
    trace = T,
    )
