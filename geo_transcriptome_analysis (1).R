# ABRINDO PACOTES NECESSÁRIOS

lapply(list("GEOquery", # baixar dados do GEO
            "edgeR", # análise estatística
            "tidyverse", # manipular dados, fazer graficos, etc ...
            "cowplot", # raincloud plot
            "ggvolc",  # volcano plot
            "ggdist", #extensão do ggplot2 que organiza distribuição de estimativas
            "ggtext", ##extensão do ggplot2 para formatação de textos
            "colorspace", #manipulação e operações associada a cores
            "PupillometryR",#manipulação de dados para aplicar em plots 
            "gghalves"), #plota gráficos que possuem mais de um elemento, ex: jitter + boxplot + etc
       require,
       character.only = T)

install.packages("ggtext")

# FAZENDO DOWNLOAD DOS DADOS DO EXPERIMENTO

raw_data <- getGEO("GSE105261", GSEMatrix = T)

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

# obtendo características das probes
feature_data <- raw_data@featureData@data

# filtrando IDs e nomes dos genes
feature_data <- feature_data[, c("ID", "ILMN_Gene")] 
# atribuindo NA aos valores vazios
feature_data[feature_data == ""] <- NA


### SUBSTITUINDO PROBES IDS PELOS NOMES DOS GENES

# obtendo dados de expressão
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

  
### SELECIONANDO PACIENTES DE INTERESSE

processed_expression_data <- expression_data %>%
  # atribuindo genes aos nomes das linhas
  `row.names<-`(.$gene) %>%
  # excluindo coluna gene
  dplyr::select(- gene) %>%
  # filtrando pacientes selecionados nos metadados
  dplyr::select(metadata$samples)
  
# checando se todas amostras dos metadados são as mesmas dos dados de expressão
all(metadata$samples %in% colnames(processed_expression_data))

# checando se todas amostras dos metadados estão na mesma ordem dos dados de expressão
all(metadata$samples == colnames(processed_expression_data))


### ANÁLISE DE EXPRESSÃO DIFERENCIAL

# criando um objeto dgeList para análise com o pacote edgeR
dge <- DGEList(processed_expression_data)

# adicionando os grupos ao dgeList como fator
dge$samples$group <- factor(metadata$tissue)

# checando diferenças nos grupos selecionados para analise (gráfico PCA)
plotMDS(dge, col = metadata$color, pch=16)

# converter numero de counts em Counts Per Million (cpm)
mycpm <- cpm(processed_expression_data)

# filtrando amostras que tenham cpm muito baixo (cpm < 0.5)
  # thresh <- mycpm > .5

# convertendo cpm em log2 cpm
logcounts <- cpm(mycpm,log=TRUE)

# checando a distribuição dos dados de contagem
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)


### NORMALIZAR DADOS

# normaliza os dados de contagem do objeto dgelist (TMM)
dge_normalizado <- calcNormFactors(dge)

# plota a distribuição da expressão genica em um paciente
plotMD(dge_normalizado,column = 7)


### EXPRESSÃO DIFERENCIAL

# Fazer matriz dummy para os grupos de estudo
design <- model.matrix(~ 0 + metadata$tissue)
# Renomeando colunas para o nome dos grupos de estudo
colnames(design) <- levels(metadata$tissue) 

# Criando um objeto Voom (transforma os dados em logCPM levando em consideração 
# a variancia media da amostra) para ser analisado pela função lmFit
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


### MANIPULANDO DADOS DE EXPRESSÃO

head(rownames_to_column(processed_expression_data, var = "x"))

graph_table <- processed_expression_data %>%
  filter(rownames(.) %in% c("COL3A1", "FAS", "TMEM119", "SEC24A", "ALDH1A3")) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "samples") %>%
  tidyr::pivot_longer(names_to = 'genes', 
                      values_to = "expression", 
                      - samples) %>%
  left_join(metadata, by = "samples") %>%
  select(-color)


expression_data_final <- expression_data %>%
  # reformatando dados para três colunas (gene [inaltarada], amostras e expressão)
  tidyr::pivot_longer(names_to = 'samples', 
                      values_to = "expression", 
                      - gene) %>%
  # juntando dados de expressão com metadados
  dplyr::left_join(metadata, 
                   by = "samples")
  
### PLOTAGEM RAINCLOUDPLOT

 #parâmetros de tamanho
  w = 8
  h = 4

raincloud2 <-
  
pal <- c( "#A034F0", "#159090") # paleta de cores
  
baseplot <-  
  ggplot(graph_table, aes(x = genes, y = expression, fill = tissue)) +
  geom_flat_violin(
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
  
  
baseplot_axis <- 
  baseplot +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Gene names",
    y = "RNA counts",
    title = "Transcriptome analysis of primary and metastasis ccRCC",
    subtitle = "Distribution of the most significantly differentially expressed genes among primary and metastatic ccRCC tumors",
    caption = "Nam HY,  *et al.* (2019) *Mol Cancer Res* DOI: 10.1158/1541-7786<br>Visualization: Vital PS and Siqueira AP"
  ) 


 baseplot_axis +
  theme_minimal(base_family = "Zilla Slab", base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), # exclude grids
    panel.grid.major = element_blank(), # exclude grids
    axis.ticks = element_blank(),
    axis.text.x = element_text(family = "Roboto Mono"),
    axis.text.y = element_text(
      color = "gray40", 
      size = 12),
    axis.title.x = element_text(
      margin = margin(t = 10),
      size = 16,
      color = "gray40"),
    axis.title.y = element_text(
      margin = margin(t = 10),
      size = 16,
      color = "gray40"),
    plot.title = element_markdown(face = "bold", size = 21),
    plot.subtitle = element_text(
      color = "grey40", hjust = 0,
      margin = margin(0, 0, 20, 0)
    ) ,
    plot.title.position = "plot",
    plot.caption = element_markdown(
      color = "grey40", lineheight = 1.2,
      margin = margin(20, 0, 0, 0)),
    plot.margin = margin(15, 15, 10, 15)
  )

ggsave('raincloudplot.png', width = 20, height = 10)

