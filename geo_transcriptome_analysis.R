### ABRINDO PACOTES NECESSÁRIOS

lapply(list("tidyverse", # manipular dados, fazer graficos, etc ...
            "ggvolc",  # volcano plot
            require,
            character.only = T))




### FAZENDO DOWNLOAD DOS DADOS DO EXPERIMENTO

raw_data <- GEOquery::getGEO("GSE105261", GSEMatrix = T)

raw_data <- raw_data[[1]]



### RECUPERANDO METADADOS

# obtendo metadados do arquivo GEO
raw_metadata <- raw_data@phenoData@data

metadata <- raw_metadata %>% 
  # filtrar metadados importantes
  dplyr::select(c("geo_accession", "tissue:ch1")) %>%
  # renomear colunas
  dplyr::rename(samples = 'geo_accession') %>%
  dplyr::rename(tissue = 'tissue:ch1') %>%
  # renomar observações da coluna tecido
  dplyr::mutate(tissue = dplyr::case_match(
    tissue,
    "Benign tissue" ~ "benign",
    "Normal renal tissue" ~ "normal",
    "Metastatic ccRCC tissue" ~ "metastatic",
    "Primary ccRCC tissue" ~ "primary")) %>%
  # filtrar tecidos normais e tumores primários
  dplyr::filter(tissue %in% c("metastatic", "primary")) %>%
  # transformando coluna tecido para fator
  dplyr::mutate(tissue = base::factor(.$tissue)) %>%
  # adicionar identificadores para a coluna tecido
  dplyr::mutate(color = dplyr::case_match(
    tissue,
    "metastatic" ~ "blue",
    "primary" ~ "green"))

# remove raw_metadata object
base::rm("raw_metadata")


### RECUPERANDOS PROBES IDS

# obtendo características das probes
feature_data <- raw_data@featureData@data

# filtrando IDs e nomes dos genes
feature_data <- feature_data[, c("ID", "ILMN_Gene")] 
# atribuindo NA aos valores vazios
feature_data[feature_data == ""] <- NA



### SUBSTITUINDO PROBES IDS PELOS NOMES DOS GENES

# obtendo dados de expressão
expression_data <- Biobase::exprs(raw_data) %>%
  # transforma para dataframe
  base::as.data.frame() %>%
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
  stats::na.omit()
  


### SELECIONANDO PACIENTES DE INTERESSE

filtered_expression_data <- expression_data %>%
  # atribuindo genes aos nomes das linhas
  base::`row.names<-`(.$gene) %>%
  # excluindo coluna gene
  dplyr::select(- gene) %>%
  # filtrando pacientes selecionados nos metadados
  dplyr::select(metadata$samples)
  
# checando se todas amostras dos metadados são as mesmas dos dados de expressão
all(metadata$samples %in% colnames(filtered_expression_data))

# checando se todas amostras dos metadados estão na mesma ordem dos dados de expressão
all(metadata$samples == colnames(filtered_expression_data))



### CHECANDO DISTRIBUIÇÃO DOS DADOS

# criando um objeto dgeList para análise com o pacote edgeR
dge_object <- edgeR::DGEList(filtered_expression_data)

# adicionando os grupos ao dgeList como fator
dge_object$samples$group <- base::factor(metadata$tissue)

# checando diferenças nos grupos selecionados para analise (gráfico PCA)
limma::plotMDS(dge_object, 
               col = metadata$color, 
               pch=16)

# converter numero de counts em Counts Per Million (cpm)
mycpm <- edgeR::cpm(filtered_expression_data)

# filtrando amostras que tenham cpm muito baixo (cpm < 0.5)
  # thresh <- mycpm > .5

# convertendo cpm em log2 cpm
logcounts <- edgeR::cpm(mycpm, log=TRUE)

# checando a distribuição dos dados de contagem
BiocGenerics::boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)

# remove objetos mycpm e logcounts
rm(list = c("mycpm", "logcounts"))



### NORMALIZAR DADOS

# normaliza os dados de contagem do objeto dgelist (TMM)
normalized_dge_object <- edgeR::calcNormFactors(dge_object)

# checa a distribuição da expressão genica normalizada em uma amostra
limma::plotMD(normalized_dge_object, 
              column = 7)



### EXPRESSÃO DIFERENCIAL

# fazer matriz dummy para os grupos de estudo
design <- stats::model.matrix(~ 0 + metadata$tissue)
# renomeando colunas para o nome dos grupos de estudo
colnames(design) <- base::levels(metadata$tissue) 


# criando um objeto Voom (transforma os dados em logCPM levando em consideração 
# a variancia media da amostra) para ser analisado pela função lmFit
voom_object <- limma::voom(normalized_dge_object,
                          design = design,
                          plot = TRUE)

# ajustar modelo linear
fit <- limma::lmFit(voom_object)

# cria uma matriz de comparação entre os grupos de estudo
contrast_matrix <- limma::makeContrasts(primary_vs_metastatic = metastatic - primary,
                                        levels = design)

# computa coeficientes e erros padrões em um modelo linear, dada uma matriz de comparação
statistical_analysis <- limma::contrasts.fit(fit = fit, 
                                 contrasts = contrast_matrix)

# análisess estatísticas (T, F, )
statistical_analysis <- limma::eBayes(statistical_analysis)

# calcula quais genes estão diferencialmente expressos, dados os contrastes
summa.fit <- limma::decideTests(statistical_analysis)

# retorna quantos genes estão diferencialmente expressos (up, down ou não significativo)
base::summary(statistical_analysis)

pvalues <- statistical_analysis$p.value

# arredonda valores da tabela para 3 digitos pós virgula
options(digits=3)

# gera tabela com valores de logFC, expressão média, t, pvalor, pvalor ajustado
top_genes_table <- limma::topTable(statistical_analysis)

# plota um volcano plot dos top 10 genes diferencialmente expressos
top_genes_table %>%
  dplyr::rename(pvalue = P.Value) %>%
  dplyr::rename(log2FoldChange = logFC) %>%
  ggvolc::ggvolc(add_seg = T )
       


### MANIPULANDO DADOS DE EXPRESSÃO


graph_table <- filtered_expression_data %>%
  # filtra os dados de expressão dos top 10 genes
  dplyr::filter(rownames(.) %in% row.names(top_genes_table)) %>%
  # transpoe a tabela
  base::t() %>%
  # transforma para data frame
  base::as.data.frame() %>%
  # transforma o nome das linhas e os transforma para uma coluna chamada samples
  tibble::rownames_to_column(var = "samples") %>%
  # transforma as variaveis samples, genes e expression em colunas mantendo samples como padrão
  tidyr::pivot_longer(names_to = 'genes', 
                      values_to = "expression", 
                      - samples) %>%
  # acrescenta os dados clínicos dos pacientes
  dplyr::left_join(metadata, by = "samples") %>%
  # exclui a coluna color
  dplyr::select(-color)



### PLOTANDO O RAINCLOUD PLOT
  
# paleta de cores
pal <- c("#A034F0", "#159090") # paleta de cores

# plot base
raincloud_plot <-
  ggplot2::ggplot(graph_table, aes(x = genes, y = expression, fill = tissue)) +
  PupillometryR::geom_flat_violin(
    position = position_nudge(x = .1, y = 0),
    adjust = 1.5,
    trim = FALSE,
    alpha = .5,
    colour = NA
  ) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = .5,
    width = .1,
    colour = "grey" # contorno do boxplot
  ) +
  geom_point(
    aes(colour = tissue),
    position = position_jitter(width = .03, height = 0),
    alpha = 0.3,
    shape = 20
  )

# plot mais nome dos eixos
raincloud_plot <-
  raincloud_plot +
  ggplot2::scale_color_manual(values = pal) +
  ggplot2::scale_fill_manual(values = pal) +
  ggplot2::labs(
    x = "Gene names",
    y = "RNA counts",
    title = "Transcriptome analysis of primary and metastasis ccRCC",
    subtitle = "Distribution of the most significantly differentially expressed genes among primary and metastatic ccRCC tumors",
    caption = "Nam HY,  *et al.* (2019) *Mol Cancer Res* DOI: 10.1158/1541-7786<br>Visualization: Vital PS and Siqueira AP"
  )

# plot mais personalizações
raincloud_plot <-
  raincloud_plot +
  ggplot2::theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    # exclude grids
    panel.grid.major = ggplot2::element_blank(),
    # exclude grids
    axis.ticks = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(family = "Roboto Mono"),
    axis.text.y = ggplot2::element_text(color = "gray40",
                                        size = 12),
    axis.title.x = ggplot2::element_text(
      margin = margin(t = 10),
      size = 16,
      color = "gray40"
    ),
    axis.title.y = ggplot2::element_text(
      margin = margin(t = 10),
      size = 16,
      color = "gray40"
    ),
    plot.title = ggtext::element_markdown(face = "bold", size = 21),
    plot.subtitle = ggplot2::element_text(
      color = "grey40",
      hjust = 0,
      margin = margin(0, 0, 20, 0)
    ) ,
    plot.title.position = "plot",
    plot.caption = ggtext::element_markdown(
      color = "grey40",
      lineheight = 1.2,
      margin = margin(20, 0, 0, 0)
    ),
    plot.margin = margin(15, 15, 10, 15)
  )

# plota o raincloud plot
raincloud_plot

# salvar grafico como imagem
ggsave('raincloudplot.png', width = 20, height = 10)
