median_Z <- function(causedata){
  library(cause)
  summary(causedata)$quants[[1]][1,3]
}


annotate_snp_loc <- function(rsid,...){


  con <- RMySQL::dbConnect(RMySQL::MySQL(),
                             host="genome-mysql.cse.ucsc.edu",
                             username="genome",
                             dbname="hg19")
  variant_df <- tbl(con,"snp150")
  anno_df <- dplyr::filter(variant_df,name %in% rsid) %>%
    dplyr::select(name,chrom,pos=chromEnd,observed,func) %>% #Pull out the SNPs we need
    dplyr::collect() %>% #grab the whole thing from the database
    dplyr::mutate(chrom=str_replace(chrom,pattern = "chr([0-9]+)",replacement = "\\1")) %>%
    dplyr::filter(!str_detect(chrom,"_")) %>%
    tidyr::separate(observed,into=c("ref","alt"),sep="/",extra="drop")

  RMySQL::dbDisconnect(con)
  functional_categories <- c("nonsense","missense","ncRNA")

  anno_df <- mutate(anno_df, has_func=map_lgl(str_split(func,","), ~any(.x %in% functional_categories)))


  return(anno_df)

}


retrieve_genes <- function(getSymbol=T){

  con <- RMySQL::dbConnect(RMySQL::MySQL(),
                             host="genome-euro-mysql.soe.ucsc.edu",
                             username="genome",
                             dbname="hg19")

  nc_genelink_df <- tbl(con,"knownToRefSeq")
  nc_genelink_df
  entrez_genelink_df <- tbl(con,"knownToLocusLink") %>% rename(ucsc_name=name,entrez_gene=value)
  kt_ensembl <- tbl(con,"knownToEnsembl")
  egene <- tbl(con,"ensGene")
  es_genelink_df <- inner_join(kt_ensembl,egene,by=c("value"="name")) %>%
    dplyr::select(ucsc_name=name,ensembl_transcript=value,ensembl_gene_stable_id=name2) %>%
    dplyr::select(-ensembl_transcript) %>%
    inner_join(entrez_genelink_df)
  kgene <- tbl(con,"knownGene")

  if(getSymbol){
    ks <- tbl(con,"kgXref")
    kgene <- inner_join(kgene,ks,by=c("name"="kgID")) %>%
      dplyr::select(ucsc_name=name,
                  chrom,
                  geneSymbol,
                  gene_start=txStart,
                  gene_end=txEnd)
  }else{
    kgene <- kgene %>%
    dplyr::select(ucsc_name=name,
                  chrom,
                  gene_start=txStart,
                  gene_end=txEnd)
  }

  gene_df <- kgene %>%
    inner_join(es_genelink_df) %>%
    dplyr::select(-ucsc_name)  %>%
    collect() %>%
    dplyr::filter(chrom %in% paste("chr",1:22,sep="")) %>%
    dplyr::mutate(chrom=stringr::str_replace(chrom,pattern = "chr([0-9]+)",replacement = "\\1")) %>%
    rename(gene_name=entrez_gene) %>% distinct(gene_name,.keep_all = T)
  RMySQL::dbDisconnect(con)
  return(gene_df)
}

unique_col <- function(...,col="snp"){
  al <- list(...)
  unique(unlist(map(al,col)))
}

nearest_gene <- function(anno_df,gene_df){
  snp_chrom <- distinct(anno_df,chrom)
  gene_chrom <- distinct(gene_df,chrom)

  snpgene_chrom <- inner_join(snp_chrom,gene_chrom)

  genemap <- map_df(snpgene_chrom$chrom,function(ch){
    dplyr::inner_join(
      filter(anno_df,chrom==ch),
      filter(gene_df,chrom==ch),by="chrom")  %>%
      dplyr::mutate(gene_dist=ifelse(pos>=gene_start & pos <=gene_end,
                                     0,
                                     .Machine$integer.max)) %>%
      dplyr::mutate(gene_dist=pmin(gene_dist,pmin(abs(pos-gene_end),
                                                  abs(pos-gene_start)))) %>%
      dplyr::group_by(name) %>%
      dplyr::filter(gene_dist==min(gene_dist)) %>%
      dplyr::distinct(gene_name,.keep_all=T) %>%  ungroup()
  })
  return(genemap)
}

retrieve_go_ids <- function(all_genes){
  ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")

  goids = biomaRt::getBM(attributes = c('ensembl_gene_id', 'go_id','name_1006','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = all_genes,
                         mart = ensembl)
}


read_eqtl_f <- function(file_name){
tiss <- fs::path_ext_remove(file_name)
read_feather(file_name,columns = c("chrom","ensembl_gene_stable_id","pos")) %>%
  mutate(tiss=tiss)
}



subset_eqtl <- function(anno_df,t_eqtl_df){
  # t_tiss_df <- filter(all_tiss_df,size==min(size))
  dplyr::distinct(anno_df,chrom,pos,name) %>%
    dplyr::inner_join(t_eqtl_df,by=c("chrom","pos"))

  return(all_eqtl)
}


finalize_snp_gene <- function(anno_df,eqtl_df,goids){
  gene_dist_df <- dplyr::select(anno_df,name,gene=ensembl_gene_stable_id)
  s_eqtl_df <- dplyr::select(eqtl_df,gene=ensembl_gene_stable_id,name,func) %>%
    dplyr::group_by(gene,name) %>%
    dplyr::summarise(ct=n_distinct(func)) %>%
    dplyr::ungroup()

  tot_df <- bind_rows(s_eqtl_df,gene_dist_df)
  dplyr::group_by(tot_df,gene,name) %>%
    dplyr::summarise(ct=n()) %>%
    dplyr::arrange(name,desc(ct)) %>% dplyr::ungroup() %>%
    dplyr::inner_join(goids,by=c("gene"="ensembl_gene_id")) %>%
    dplyr::mutate(is_pathway=T) %>%
    dplyr::filter(go_id!="")
}

split_slice <- function(.data,num_slices=10){

  spliti <- parallel::splitIndices(nrow(.data),num_slices)
  map(spliti,~slice(.data,.x))
}


three_models <- function(prob_Z1_conf,is_pathway,is_shared){
  dplyr::bind_rows(dplyr::mutate(broom::tidy(lm(is_pathway~prob_Z1_conf)),method="lm"),
            dplyr::mutate(broom::tidy(glm(is_pathway~prob_Z1_conf,family="binomial")),method="logit"),
            dplyr::mutate(broom::tidy(fisher.test(is_shared,is_pathway)),method="fe",term="prob_Z1_conf") %>%
              dplyr::select(-conf.low,-conf.high,-alternative)
  )
}

enrich_fun <- function(cause_df,mz,snp_gene_df,file,go_id_df=NULL){
  if(!is.null(go_id_df)){
    snp_gene_df <- dplyr::semi_join(snp_gene_df,go_id_df)
  }

  sub_z <- set_names(cause_df$prob_Z1_conf,cause_df$snp)
  td <-   as_tibble(cause_df) %>%
    dplyr::select(snp,prob_Z1_conf) %>%
    dplyr::inner_join(snp_gene_df,by=c("snp"="name")) %>%
    dplyr::select(snp,gene,prob_Z1_conf,go_id,is_pathway) %>%
    dplyr::mutate(snp=factor(snp)) %>%
    dplyr::distinct(snp,prob_Z1_conf,go_id,is_pathway)
  if(nrow(td)>0){

    ret_df <- td %>% dplyr::group_by(go_id) %>%
      tidyr::complete(snp,go_id,fill=list(is_pathway=FALSE)) %>%
      dplyr::mutate(prob_Z1_conf=sub_z[as.character(snp)],is_shared=prob_Z1_conf>mz) %>%
      dplyr::distinct() %>%
      dplyr::do(three_models(.$prob_Z1_conf,.$is_pathway,.$is_shared)) %>%
      dplyr::ungroup() %>% dplyr::group_by(method) %>%
      dplyr::mutate(p_BH=p.adjust(p.value,method="BH"),p_Bon=p.adjust(p.value,method = "bonferroni")) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(p.value) %>% dplyr::mutate(file=file)
    return(ret_df)

  }
  else{
    return(NULL)
  }
}



