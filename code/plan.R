
c_dir <- fs::path_expand("~/Desktop/cause_dir/")

eqtl_dir <- fs::path_expand("/run/media/nwknoblauch/Data/eQTL_SumStats/GTEx_Analysis_v7_eQTL/")
all_files <- fs::path_file(fs::dir_ls(c_dir,glob="*cause.RDS"))
eqtl_files <- fs::path_file(fs::dir_ls(eqtl_dir,glob="*fst"))
plan <- drake_plan(
  causedata = target(readRDS(fs::path(c_dir,all_files[x])),transform = map(x=!!(seq_along(all_files)))),
  medz = target( median_Z(causedata),transform=map(causedata)),
  cause_df = target( causedata$data,transform=map(causedata)),
  all_snps = target(unique_col(cause_df,col="snp"),transform =combine(cause_df)),
  snp_df = target(annotate_snp_loc(all_snps),hpc=FALSE),
  gene_df = target(retrieve_genes(),hpc=FALSE),
  anno_df = nearest_gene(snp_df,gene_df = gene_df),
  t_eqtl_df = target(read_eqtl_f(fs::path(eqtl_dir,eqtl_files[i])),transform=map(i=!!seq_along(eqtl_files))),
  s_eqtl_df = target(subset_eqtl(anno_df,t_eqtl_df = t_eqtl_df),transform=map(t_eqtl_df)),
  eqtl_df = target(bind_rows(s_eqtl_df),transform=combine(s_eqtl_df)),
  genelist = unique_col(eqtl_df,anno_df,col="ensembl_gene_stable_id"),
  godf = retrieve_go_ids(genelist),
  snp_gene_df =finalize_snp_gene(anno_df,eqtl_df,godf),
  all_go_terms = tibble::tibble(go_id=unique(snp_gene_df$go_id)),
  split_go_df = split_slice(all_go_terms,10),
  res_df = target(enrich_fun(cause_df,medz,snp_gene_df,file=all_files[x]), transform = map(cause_df,medz,x=!!(seq_along(all_files)))),

  max_expand=4L,
  trace=T,
)
