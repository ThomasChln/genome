render_datatable = function(df) {
  col_defs = c('details', 'study', 'initial.sample.size', 'replication.sample.size') %>% {
      list(list(width = '400px', targets = match(., names(df)) - 1))
    }

  buttons = c('copy', 'csv', 'excel', 'pdf', 'print') %>%
    list(text = 'Export', extend = 'collection', buttons = .) %>%
    list('pageLength', 'colvis', .)

  DT::datatable(df, escape = FALSE, rownames = FALSE, filter = 'top',
    extensions = 'Buttons',
    options = list(pageLength = 25, dom = 'Bfrtip', buttons = buttons,
      autoWidth = TRUE, columnDefs = col_defs))
}

tidy_df = function(obj) {
  cnames = c('snps', 'genotype', 'magnitude', 'repute', 'summary', 'details',
    'snpedia_url', 'link', 'date', 'study', 'disease.trait',
    'risk.allele.frequency', 'or.or.beta', 'initial.sample.size',
    'replication.sample.size', 'risk_allele', 'pvalue_mlog', 'x95..ci..text.')

  obj %>% setNames(tolower(names(.))) %>% `[`(cnames) %>%

    # fix urls
    `$<-`('snpedia_url', .$snpedia_url %>% 
        paste0('https://www.snpedia.com/index.php/', .) %>%
        paste0('<a href="', ., '">', ., '</a>')) %>%
    `$<-`('link', .$link %>% paste0('https://', .) %>%
        paste0('<a href="', ., '">', ., '</a>')) %>%

    # format columns
    `$<-`('pvalue_mlog', round(.$pvalue_mlog, 3)) %>%
    `$<-`('magnitude', .$magnitude %>% as.character %>% as.numeric) %>%
    `$<-`('genotype', .$genotype %>% factor) %>%

    # rename columns
    `names<-`(gsub('link', 'EBI study', names(.))) %>%
    `names<-`(gsub('snps', 'marker', names(.))) %>%

    `[`(order(.$magnitude, decreasing = TRUE), )
}

# risk allele orientation depending of author ? neither always dbsnp nor (+)
search_risk_allele = function(snp) {
  if (snp['OR.or.BETA'] %>% as.numeric < 1) {
    snp['risk_allele'] %<>% paste0('[^', ., ']')
  }
  grepl(snp['risk_allele'], snp['genotype'])
}

null_omit = function(obj) obj %>% `[`(!sapply(., is.null))

check_strand = function(allele, strand) {
  allele %>% switch(strand, plus = .,
      minus = switch(., A = 'T', T = 'A', C = 'G', G = 'C'))
}

format_genotype = function(geno, strand) {
  sapply(1:2, function(i) {
      substr(geno, i, i) %>% check_strand(strand)
    }) %>% sort %>% {
      paste0('(', .[1], ';', .[2], ')')
    }
}

snpedia_url = function(snp) {
  if (is.na(snp['Orientation'])) return(NA)
  format_genotype(snp['genotype'], snp['Orientation']) %>% 
    paste0(snp['SNPS'], .)
}

get_genotype_details = function(summary) {
  gsub('.*}}\n?', '', summary) %>% gsub('\n+', ' ', .)
}

format_genotype_answer = function(summaries) {
  sapply(summaries, SNPediaR::extractTags,
    c('rsid', 'magnitude', 'repute', 'summary')) %>% t %>%
    cbind.data.frame(snpedia_url = names(summaries),
      details = sapply(summaries, get_genotype_details))
}

snpedia_snp_tags = function(df) {
  df$SNPS %>% SNPediaR::getPages() %>% null_omit %>%
    sapply(SNPediaR::extractSnpTags, c('rsid', 'Orientation')) %>% t
}

snpedia_genotype_tags = function(df) {
  df %>% apply(1, snpedia_url) %>% na.omit %>% SNPediaR::getPages() %>%
    null_omit %>%  format_genotype_answer
}

qc_ebi = function(df, pval_mlog_threshold) {
  subset(df, !is.na(OR.or.BETA) & OR.or.BETA != 0 &
      PVALUE_MLOG != Inf & PVALUE_MLOG > pval_mlog_threshold) %>%

    cbind(risk_allele = gsub('.*-', '', .$STRONGEST.SNP.RISK.ALLELE)) %>%
    subset(risk_allele %in% c('A', 'C', 'T', 'G')) %>%

    `[`(order(.$PVALUE_MLOG, decreasing = TRUE), ) %>%
    subset(!duplicated(SNPS))
}

ebi_risk_snps = function(ebi_path = 'ebi_gwas_catalog.tsv',
  pval_mlog_threshold = 10) {

  if (!file.exists(ebi_path)) {
    'https://www.ebi.ac.uk/gwas/api/search/downloads/full' %>%
      download.file(ebi_path)
  }
  read.csv(ebi_path, sep = '\t') %>% qc_ebi(pval_mlog_threshold)
}

get_snpedia = function(df) {
  df %>% merge(snpedia_snp_tags(.), by.x = 'SNP_ID_CURRENT', by.y = 'rsid',
      all.x = TRUE) %>%
    merge(snpedia_genotype_tags(.), ., by.x = 'rsid', by.y = 'SNP_ID_CURRENT',
      all.y = TRUE)
}

get_ebi = function(df) {
  df %>% merge(ebi_risk_snps(), ., by.x = 'SNPS', by.y = '# rsid') %>%
    subset(apply(., 1, search_risk_allele))
}

get_genome = function(path) {
  path %>% ifelse(grepl('.zip$', .), unzip(., exdir = '/tmp'),  .) %>%
    data.table::fread(data.table = FALSE) %>%
    subset(nchar(genotype) == 2 & !grepl('D|I|-', genotype))
}

get_summary = function(genome) genome %>% get_genome %>% get_ebi %>% get_snpedia

load_or_compute = function(.) {
  hash = paste0(., '_sha256_', digest::digest(., 'sha256'), '.rds')
  if (file.exists(hash)) {
    get(load(hash))
  } else {
    get_summary(.) %T>% save(file = hash)
  }
}
