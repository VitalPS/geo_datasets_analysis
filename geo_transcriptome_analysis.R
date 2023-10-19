# ABRINDO PACOTES NECESSÁRIOS

lapply(list("GEOquery", # baixar dados do GEO
            "DESeq2", # análise estatística
            "edgeR", # análise estatística
            "tidyverse", # manipular dados)
            "ggplot2", # fazer gráficos
            "ggvolc"
            ), 
       require,
       character.only = T)



# FAZENDO DOWNLOAD DOS DADOS DO EXPERIMENTO

raw_data <- getGEO("GSE105261", GSEMatrix = T)

raw_data <- raw_data[[1]]


geo <- getGEO("GSE105261", GSEMatrix = T)

geo <- raw_data[[1]]




### RECUPERANDO METADADOS

raw_metadata <- raw_data@phenoData@data

metadata <- raw_metadata %>% 
  # filtrar metadados importantes
  dplyr::select(c("geo_accession", "tissue:ch1")) %>%
  # renomear colunas
  dplyr::rename(samples = 'geo_accession') %>%
  dplyr::rename(tissue = 'tissue:ch1') %>%
  # renomar observações da coluna tecido
  mutate(tissue = case_match(
    tissue,
    "Benign tissue" ~ "benign",
    "Normal renal tissue" ~ "normal",
    "Metastatic ccRCC tissue" ~ "metastatic",
    "Primary ccRCC tissue" ~ "primary")) %>%
  # filtrar tecidos normais e tumores primários
  dplyr::filter(tissue %in% c("metastatic", "primary")) %>%
  # transformando coluna tecido para fator
  mutate(tissue = factor(.$tissue)) %>%
  # adicionar identificadores para a coluna tecido
  mutate(color = case_match(
    tissue,
    "metastatic" ~ "blue",
    "primary" ~ "green"))




### RECUPERANDOS PROBES IDS

feature_data <- raw_data@featureData@data

# filtrando IDs e nomes dos genes
feature_data <- feature_data[, c("ID", "ILMN_Gene")] 
# atribuindo NA aos valores vazios
feature_data[feature_data == ""] <- NA



### SUBSTITUINDO PROBES IDS PELOS NOMES DOS GENES

# Recuperando dados de expressão
raw_expression_data <- as.data.frame(exprs(raw_data))

expression_data <- raw_expression_data %>%
  # cria uma coluna com os simbolos
  tibble::rownames_to_column(var = "ID") %>% 
  # juntando probes ids com dados de expressão
  dplyr::inner_join(feature_data, by = "ID") %>% 
  # renomeando colunas
  dplyr::rename(gene = ILMN_Gene) %>%
  # removendo colunas não importantes
  dplyr::select(-ID) %>%
  # removendo itens duplicados
  dplyr::filter(!duplicated(.$gene)) %>%
  # removendo os NAs
  na.omit()

  


### NORMALIZAR E FAZER ANÁLISE ESTATÍSTICA

processed_expression_data <- expression_data %>%
  `row.names<-`(.$gene) %>%
  dplyr::select(- gene) %>%
  dplyr::select(metadata$samples)
  


# checando se todas amostras dos metadados são as mesmas dos dados de expressão
all(metadata$samples %in% colnames(processed_expression_data))


# checando se todas amostras dos metadados estão na mesma ordem dos dados de expressão
all(metadata$samples == colnames(processed_expression_data))


### ANÁLISE ESTATÍSTICA

# criando um objeto dgeList
dge <- DGEList(processed_expression_data)

# adicionando os grupos ao dgeList como fator
dge$samples$group <- factor(metadata$tissue)

## Plotar gráfico PCA
plotMDS(dge, col = metadata$color, pch=16)




mycpm <- cpm(processed_expression_data)

thresh <- mycpm > .5

logcounts <- cpm(mycpm,log=TRUE)

boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)



### NORMALIZAR DADOS

dge_normalizado <- calcNormFactors(dge)

plotMD(dge_normalizado,column = 7)




### EXPRESSÃO DIFERENCIAL

# Fazer matriz dummy para os grupos
design <- model.matrix(~ 0 + metadata$tissue)
# Renomeando colunas
colnames(design) <- levels(metadata$tissue) 


# Criando um objeto Voom
voom_object <- voom(dge_normalizado,
                    design = design,
                    plot = TRUE)

# Modelo linear
fit <- lmFit(voom_object)

# Cria uma matriz de comparação entre diferentes grupos
cont.matrix <- makeContrasts(primary_vs_metastatic = metastatic - primary,
                              levels = design)

# 
fit.cont <- contrasts.fit(fit, cont.matrix)

fit.cont <- eBayes(fit.cont)

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

pvalues <- fit.cont$p.value

# arredonda valores da tabela para 3 digitos pós virgula
options(digits=3)

# gera tabela com valores de logFC, expressão média, t, pvalor, pvalor ajustado
final_table <- topTable(fit.cont)

final_table_processed <- final_table %>%
  rename(pvalue = P.Value) %>%
  rename(log2FoldChange = logFC)



ggvolc(final_table_processed,
       add_seg = T )
       
plot %>%
  genes_table(final_table_processed)

# Criando um dataset de expressão diferencial
# dds <- DESeqDataSetFromMatrix(
#   countData = round(processed_expression_data),
#   colData = metadata,
#   design = ~ tissue)

# dds <- estimateDispersionsGeneEst(dds)
# 
# dispersions(dds) <- mcols(dds)$dispGeneEst

# Setando controle basal
# deseq_dataset$tissue <- relevel(
#   deseq_dataset$tissue, 
#   ref = "primary"
#   )


# Executando DESeq
# dds <- DESeq()
# 
# dds_results < results(dds)



### MANIPULANDO DADOS DE EXPRESSÃO

expression_data_final <- expression_data %>%
  # reformatando dados para três colunas (gene [inaltarada], amostras e expressão)
  tidyr::pivot_longer(names_to = 'samples', 
                      values_to = "expression", 
                      - gene) %>%
  # juntando dados de expressão com metadados
  dplyr::left_join(metadata, 
                   by = "samples")
  








# 
# # criar uma lista de referência simbolo/nome do gene
# gene_simbols = data.frame(gene = unlist(mget(x = row.names(raw_expression_data),
#                 envir = illuminaHumanv4SYMBOL)))
# 
# # ignorar simbolos sem nomes encontrados (valores NA)
# gene_simbols <- na.omit(gene_simbols)
# gene_simbols[, "row_names"] <- row.names(gene_simbols)
# 
# # renomeando simbolos
# expression_data <- raw_expression_data %>%
#   dplyr::inner_join(
#     gene_simbols,
#     by = "row_names"
#   ) %>%
#   `row.names<-`(.$gene) %>%
#   dplyr::select(-c("row_names", "gene"))
# 
# duplicated(expression_data$gene)
# 
# 
# 
# # # log2 transform
# # qx <-as.numeric(
# #   quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), 
# #            na.rm = T))
# # 
# # LogC <- (qx[5] > 100) ||
# #   (qx[6] - qx[1] > 50 && qx[2] > 0)
# # 
# # if (LogC) {
# #   ex[which(ex <= 0)] <- NaN
# #   ex <- log2(ex)
# # }
# 
# 
# # box-and-whisker plot
# dev.new(width = 3 + ncol(gset) / 6, height = 5)
# 
# par(mar = c(7, 4, 2, 1))
# title <- paste (
#   "GSE105288", 
#   "/", 
#   annotation(gset), 
#   sep = ""
#   )
# 
# boxplot(
#   ex,
#   boxwex = 0.7,
#   notch = T,
#   main = title,
#   outline = FALSE,
#   las = 2
# )
# 
# dev.off()
# 
# # expression value distribution plot
# par(mar = c(4, 4, 2, 1))
# title <- paste (
#   "GSE105288", "/", 
#   annotation(gset), 
#   "value distribution", 
#   sep = "")
# 
# plotDensities(ex, main = title, legend = F)
# 
# 
# # mean-variance trend
# ex <- na.omit(ex) # eliminate rows with NAs
# 
# plotSA(
#   lmFit(ex), 
#   main = "Mean variance trend, GSE105288"
#   )
# 
# 
# 
# # # UMAP plot (multi-dimensional scaling)
# # ex <- ex[!duplicated(ex),]  # remove duplicates
# # ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
# # plot(
# #   ump$layout,
# #   main = "UMAP plot, nbrs=15",
# #   xlab = "",
# #   ylab = "",
# #   pch = 20,
# #   cex = 1.5
# # )
# # library("maptools")  # point labels without overlaps
# # pointLabel(
# #   ump$layout,
# #   labels = rownames(ump$layout),
# #   method = "SANN",
# #   cex = 0.6
# # )
