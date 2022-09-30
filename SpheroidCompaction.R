library(tidyverse)
library(ggtext)
library(patchwork)
library(ggplot2)
library(dplyr)


# LOAD DATASETS ###############################################################
Path <- 'G:/2021-zeiss 900/20220209/split max projection/Mergelayer/only4/'
data <- read.csv(file = 'G:/2021-zeiss 900/20220209/split max projection/Mergelayer/only4/Stiff-RGDBIGDATA.csv')
data$Channel <- as.factor(data$Channel)
data$file <- as.factor(data$file)
# levels(data$Channel) <- c("BF", "Necrotic", "Proliferating", "DAPI")
levels(data$Channel) <- c("Necrotic", "Proliferating", "DAPI")
###############################################################################
colorset <- c('BF'='grey', 'Necrotic'='#F8766D', 'Proliferating'='#00BA38', 'DAPI'='#619CFF')

data <- data[data$Channel != 'BF',]

# data <- data[data$file == "soft",]
# filter(data, (file == factor("stiff") & Channel == 3))
ggplot(data , aes(x= x, y= y, color= factor(Channel))) + geom_point() +
  # scale_color_viridis(discrete = TRUE, option = "D")+
  # scale_fill_viridis(discrete = TRUE) +
  theme_minimal()

plot_lst <- vector("list", length = length(levels(data$file)))
for (filename in levels(data$file)) {

  data_filename <- data%>% filter(file == filename)
  f_filename<- ggplot(data_filename , aes(x= x, y= y,color= factor(Channel))) +
    geom_point(alpha = 0.3) +
    scale_fill_manual(values = colorset) +
    theme_minimal()
  print(f_filename)
  plot_lst[[filename]] <- f_filename
  ggsave(plot = f_filename, path = Path, filename = filename, width = 8, height = 6, device = "png")
}

plot_grid(plotlist = plot_lst, width = 8, height = 6,
          ncol =3, align = "vh")

# DAPI_color <- "#BEBEBE"
# KI67_color <- "#0000FF"
# HIFalpha1_color <- "#FF0000"

######################################################################################
# StackBarPlots ######################################################################
# data <- read.csv(file = 'G:/2021-zeiss 900/20220209/split max projection/mergg/Sample6BIGDATA.csv')
# data$Channel <- as.factor(data$Channel)
# data$file <- as.factor(data$file)
# levels(data$Channel) <- c("BF", "Necrotic", "Proliferating", "DAPI")
# colorset = c('BF'='grey','Necrotic'='#F8766D','Proliferating'='#00BA38','DAPI'='#619CFF')
levels(data$Channel) <- c("Necrotic", "Proliferating", "DAPI")
colorset = c('Necrotic'='#F8766D','Proliferating'='#00BA38','DAPI'='#619CFF')

data <- data[data$Channel != 'BF',]
# data <- data[data$File != 'Sample2',]

summ_total <- data.frame()

for (filename in levels(data$file)) {
  data_filename <- data%>% filter(file == filename)
  summ<- as.data.frame(summary(data_filename$Channel))%>%
    mutate(percent = 100 * summary(data_filename$Channel)/sum(summary(data_filename$Channel)),
           File = filename, Channel = levels(data_filename$Channel))
  summ_total <- as.data.frame(rbind(summ_total,summ))
}
# summ_total$Channel <- as.factor(summ_total$Channel, levels=c("BF", "Necrotic", "Proliferating", "DAPI"))

levels(data_filename$Channel)

summ_total

colnames(summ_total)

sbp<- ggplot(summ_total, aes(fill= Channel , y= percent, x=File)) +
  geom_bar(position="stack", stat="identity", alpha = 0.7) +
  theme_classic()
sbp + scale_fill_manual(values = colorset) +
  labs(x = NULL )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########################################################################################
head(data)

cells_count <-data %>%
  group_by(Channel) %>%
  filter(file == "soft") %>%
  tally()

print(cells_count[1, 2])

cells_count <- data %>%
  count(file)
print(cells_count)

necrotic_n <- data %>%
  filter( file == "soft" & Channel == "Necrotic") %>%
  pull(1)
print(length(necrotic_n))


##########################################################################################

rlang::last_error()

diarrhea_n <- disease_count %>%
  filter(disease_stat == "DiarrhealControl") %>%
  pull(n)

case_n <- disease_count %>%
  filter(disease_stat == "Case") %>%
  pull(n)


kt <- kruskal.test(invsimpson ~ disease_stat, data=metadata_alpha)

if(kt$p.value < 0.05){
  pt <- pairwise.wilcox.test(metadata_alpha$invsimpson,
                             g=metadata_alpha$disease_stat,
                             p.adjust.method = "BH")
}

strip_chart <- metadata_alpha %>%
  ggplot(aes(x=disease_stat, y=invsimpson, fill=disease_stat)) +
  stat_summary(fun = median, show.legend=FALSE, geom="crossbar") +
  geom_jitter(show.legend=FALSE, width=0.25, shape=21, color="black") +
  labs(x=NULL,
       y="Inverse Simpson Index") +
  scale_x_discrete(breaks=c("NonDiarrhealControl","DiarrhealControl","Case"),
                   labels=c(glue("Healthy<br>(N={healthy_n})"),
                            glue("Diarrhea and<br>*C.difficile* negative<br>\\
                                 (N={diarrhea_n})"),
                            glue("Diarrhea and<br>*C.difficile* positive<br>\\
                                 (N={case_n})"))
  ) +
  scale_fill_manual(name=NULL,
                    breaks=c("NonDiarrhealControl","DiarrhealControl","Case"),
                    labels=c("Healthy",
                             "Diarrhea and<br>*C.difficile* negative",
                             "Diarrhea and<br>*C.difficile* positive"),
                    values=c(healthy_color, diarrhea_color, case_color)) +
  theme_classic() +
  theme(axis.text.x = element_markdown(size=6)) +
  geom_line(data=tibble(x=c(2, 3), y=c(23, 23)),
            aes(x=x, y=y),
            inherit.aes=FALSE) +
  geom_line(data=tibble(x=c(1, 2.5), y=c(33, 33)),
            aes(x=x, y=y),
            inherit.aes=FALSE) +
  geom_text(data=tibble(x=2.5, y=24),
            aes(x=x, y=y, label="n.s."),
            inherit.aes=FALSE) +
  geom_text(data=tibble(x=1.75, y=33.5),
            aes(x=x, y=y, label="*"), size=8,
            inherit.aes=FALSE)



disease_invsimpson <- inner_join(metadata, alpha_diversity,
                                 by=c('sample_id'='group')) %>%
  select(disease_stat, invsimpson)


get_roc_data <- function(negative, positive){

  disease_invsimpson %>%
    filter(disease_stat == negative | disease_stat == positive) %>%
    mutate(disease_stat = recode(disease_stat,
                                 "{negative}" := FALSE,
                                 "{positive}" := TRUE)) %>%
    mutate(sens_spec = map(invsimpson, get_sens_spec, .)) %>%
    unnest(sens_spec) %>%
    mutate(comparison = glue("{negative}_{positive}"))

}

get_sens_spec <- function(x, data){

  predicted <- x > data$invsimpson

  tp <- sum(predicted & data$disease_stat)
  tn <- sum(!predicted & !data$disease_stat)
  fp <- sum(predicted & !data$disease_stat)
  fn <- sum(!predicted & data$disease_stat)

  specificity <- tn / (tn + fp)
  sensitivity <- tp / (tp + fn)


  tibble("sensitivity" = sensitivity, "specificity"=specificity)
}

roc_data <- bind_rows(
  get_roc_data("NonDiarrhealControl", "DiarrhealControl"),
  get_roc_data("NonDiarrhealControl", "Case"),
  get_roc_data("DiarrhealControl", "Case")
) %>%
  arrange(invsimpson)


pretty_names <- c("NonDiarrhealControl_Case" = "<strong style='color:#BEBEBE;'>Healthy</strong> vs. <br><strong style='color:#FF0000;'>*C. difficile* positive</strong>",
                  "NonDiarrhealControl_DiarrhealControl" = "<strong style='color:#BEBEBE;'>Healthy</strong> vs.<br><strong style='color:#0000FF;'>*C. difficile* negative</strong>",
                  "DiarrhealControl_Case" = "*C. difficile*<br><strong style='color:#0000FF;'>negative</strong> vs. <br><strong style='color:#FF0000;'>positive</strong>")

get_roc_curve <- function(test){

  roc_data %>%
    filter(comparison == test) %>%
    ggplot(aes(x=1-specificity, sensitivity, group=comparison)) +
    geom_abline(slope=1, intercept=0, color="gray") +
    geom_line(size = 1, linejoin="round", color="black") +
    labs(x="1-Specificity",
         y="Sensitivity") +
    theme_classic() +
    theme(
      legend.text = element_markdown(),
      legend.key.height = unit(25, "pt"),
      legend.position = c(0.8, 0.2)
    ) +
    geom_richtext(data=tibble(x=0.75, y=0.15, label=pretty_names[test]),
                  aes(x=x, y=y, label=label), inherit.aes = FALSE,
                  fill=NA, label.color=NA)
}

ndc_c <- get_roc_curve("NonDiarrhealControl_Case")
dc_c <- get_roc_curve("DiarrhealControl_Case")
ndc_dc <- get_roc_curve("NonDiarrhealControl_DiarrhealControl")

plot_grid(strip_chart, ndc_dc, ndc_c, dc_c,
          nrow=2, ncol=2, rel_widths = c(1, 1),
          labels=c("A", "B", "C", "D"))


ggsave("schubert_fig_1.tiff", width=6.5, height=6)