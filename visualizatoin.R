aligned_dir <- '/length/'
strain <- c('CEA17', 'CEA17_dox', 'dclA_pksp_IRT', 'dclB_pksp_IRT')

require(ggplot2)
require(ggseqlogo)
require(data.table)
require(tidyverse)

first_nt <- list()
for (s in strain) {
  tmp2 <- read_tsv(paste0(aligned_dir,s, '_merge.fa_1st_nt.txt'),
                   col_names = F) 
  
  # change T to U
  tmp2$X1 <- str_replace_all(tmp2$X1, 'T', 'U')
  first_nt[s] <- tmp2
}
# Create custom colour scheme
cs1 = make_col_scheme(chars=c('A', 'U', 'C', 'G'), 
                      cols=c( "#E69F00", "#0072B2", "#D55E00", "#56B4E9"))

ggseqlogo(first_nt, 
          ncol = 3, 
          stack_width = 0.9,
          method = 'prob',
          col_scheme=cs1)
ggsave('1st_nt.pdf')


# ==== length ====
length_merge <- data.frame()
for (s in strain) {
  tmp <- fread(paste0(aligned_dir, s, '_merge_align.fa_length.txt')) %>% 
    mutate(strain = s) 
  
  length_merge <- rbind(tmp, length_merge)
}

length_merge %>% group_by(V1, strain) %>% 
  count()

ggplot(length_merge,
       aes(V1))+
  geom_histogram()+
  facet_wrap(~strain, ncol = 3)+
  scale_x_continuous(limits = c(10, 30))+
  scale_y_continuous(labels = comma)+
  xlab('Read Length')+
  ylab('Length Distribution')+
  theme_base()+
  theme(
    axis.title.y = element_text(margin = margin(r = 20), size = 18),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 16),
    #panel.border = element_blank(),  # Remove panel border
    panel.background = element_blank(),
    plot.background = element_blank())
ggsave('length_distribution.pdf')
