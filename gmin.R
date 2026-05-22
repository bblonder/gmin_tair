library(plantecophys)
library(readxl)
library(ggplot2)
library(dplyr)
library(ggpubr)

# make output dir
if (!file.exists('outputs'))
{
  dir.create('outputs')  
}

data_slot = read_excel('data_slot/nph17626-sup-0002-tabless1-s2.xlsx',sheet='Table_S1',skip=1) %>%
  rename(Tair_degC=`Temperature (C)`) %>%
  mutate(gmin_mol = Gmin_mmol / 1000)

data_area1 = read_excel('data_slot/LMA_Source_PCE2021.xlsx') %>% select(Species, Width_cm) %>% na.omit
data_area2 = read_excel('data_slot/LMA extra.xlsx') %>% select(Species, Width_cm) # 3 species estimated from Kew herbarium sheets

data_slot_joined = data_slot %>% 
  left_join(rbind(data_area1, data_area2) ,by='Species') %>%
  mutate(stomatal_ratio = 1) %>% # based on personal communication from martijn
  mutate(Dataset='Panama') %>%
  select(Dataset, Species, Tair_degC, Gmin_mmol, Width_cm, stomatal_ratio)



data_garen_traits = read.csv('data_garen/cond.traits.data_GarenMichaletz25_stm.csv') %>%
  mutate(stomatal_ratio = ifelse(stomatal_distribution=='hypostomatic', 1, 2)) %>%
  mutate(Dataset='Canada') %>%
  select(Dataset, Species=species_full, Tair_degC=treatment, 
         Gmin_mmol=gmin_mmol_m2s, stomatal_ratio, Width_cm = L_cm)



data_combined = rbind(data_garen_traits, data_slot_joined)


ggplot(data_combined, aes(x=Tair_degC,y=Gmin_mmol)) +
  facet_wrap(~Species+Dataset) +
  geom_point()


coef_table_combined = data_combined %>%
  group_by(Species) %>%
  #filter(Tair_degC > 30) %>%
  do(data.frame(t(coef(lm(Gmin_mmol ~ Tair_degC + I(Tair_degC^2), data = .))))) %>%
  rename(intercept=2,slope1=3, slope2=4)

# check fits
pdf(file='outputs/g_gmin_fits.pdf')
by(data_combined, data_combined$Species, function(x) {
  plot(Gmin_mmol~Tair_degC,data=x,main=x$Species[1])
  xv=20:50
  coef_table_combined_this = coef_table_combined %>% filter(Species==x$Species[1])
  lines(xv, coef_table_combined_this$intercept + coef_table_combined_this$slope1*xv+coef_table_combined_this$slope2*xv^2,col='red')
  }) 
dev.off()


process_species <- function(species_this, Tair_degC_focal)
{
  # make data frame of observed values
  Tair_this = data_combined %>% filter(Species==species_this) %>% pull(Tair_degC)
  gmin_mmol_this = data_combined %>% filter(Species==species_this) %>% pull(Gmin_mmol) 
  df_this = data.frame(Tair_degC=Tair_this,gmin_mmol=gmin_mmol_this)
  
  # get species traits
  width_cm_this = data_combined %>% filter(Species==species_this) %>% pull(Width_cm) %>% mean
  stomatal_ratio_this = data_combined %>% filter(Species==species_this) %>% pull(stomatal_ratio) %>% mean
  
  # predict gmin as a function of tair
  seq_Tair = seq(30,50,by=1)
  intercept_this = coef_table_combined %>% filter(Species==species_this) %>% pull(intercept)
  slope1_this = coef_table_combined %>% filter(Species==species_this) %>% pull(slope1)
  slope2_this = coef_table_combined %>% filter(Species==species_this) %>% pull(slope2)
  seq_gmin_mmol = (intercept_this + slope1_this*Tair_degC_focal + slope2_this*Tair_degC_focal^2)
  
  
  params = data.frame(Wind=1,
                      Wleaf = width_cm_this / 100, # convert from cm to m
                      StomatalRatio = stomatal_ratio_this,
                      LeafAbs = 0.86,
                      gs=seq_gmin_mmol / 1000, # convert from mmol to mol
                      Tair=seq_Tair,
                      VPD=RHtoVPD(RH=0.5,TdegC=seq_Tair,Pa=101))
  
  
  tleaf_with_gmin = apply(params, 1, FUN=function(x) {
    x = as.list(x)
    do.call("FindTleaf",args=x)
    })
  tleaf_without_gmin = apply(params, 1, FUN=function(x) {
    x = as.list(x)
    x$gs = 0
    do.call("FindTleaf",args=x)
    })
  
  params$tleaf_with_gmin = tleaf_with_gmin
  params$tleaf_without_gmin = tleaf_without_gmin
  params$cooling_effect_degc = params$tleaf_with_gmin - params$tleaf_without_gmin
  params$delta_t_degc = params$tleaf_with_gmin - params$Tair

  
  g1 = ggplot(params, aes(x=Tair,y=cooling_effect_degc)) +
    geom_line() +
    theme_bw()
  
  g2 = ggplot(params, aes(x=Tair,y=delta_t_degc)) +
    geom_line() +
    theme_bw()
  
  ggsave(ggarrange(g1, g2, align='hv',nrow=1,ncol=2),
         file=sprintf('outputs/result_%s.pdf',species_this),width=8,height=3.5)
  
  cat('.')
  
  return(data.frame(Species=species_this,
                    delta_t_degc = params %>% filter(Tair==Tair_degC_focal) %>% pull(delta_t_degc),
                    cooling_effect_degc = params %>% filter(Tair==Tair_degC_focal) %>% pull(cooling_effect_degc)))
}

# run all species
species_list_combined = sort(unique(data_combined$Species))
results_45 = do.call('rbind',lapply(species_list_combined, process_species, Tair_degC_focal=45)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info

g_cooling_effect = ggplot(results_45, aes(x=reorder(Species, cooling_effect_degc),y=cooling_effect_degc,fill=Dataset)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_bw() +
  ylab('Cooling effect (°C)') +
  xlab('Species') +
  geom_hline(yintercept = 0, color='black')
ggsave(g_cooling_effect,file='outputs/g_cooling_effect.pdf',width=7,height=7)

g_delta_t = ggplot(results_45, aes(x=reorder(Species, delta_t_degc),y=delta_t_degc,fill=Dataset)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_bw() +
  ylab('∆T (°C)') +
  xlab('')
ggsave(g_delta_t,file='outputs/g_delta_t.pdf',width=7,height=7)



g_pairs = ggplot(results_45, aes(x=delta_t_degc, y=cooling_effect_degc)) + 
  geom_point(color='black') +
  xlim(-1,10.5) + ylim(-2,0) +
  geom_hline(yintercept = 0,color='gray') +
  geom_vline(xintercept = 0,color='gray') +
  theme_bw() +
  xlab('∆T (°C)') +
  ylab('Cooling effect (°C)') +
  stat_smooth(method='lm',se=FALSE)
ggsave(g_pairs, file='outputs/g_pairs.pdf', width=7,height=7)


g_fig1 = ggarrange(g_cooling_effect, g_delta_t, g_pairs, labels='AUTO',nrow=2,ncol=2)
ggsave(g_fig1, file='outputs/g_fig1.pdf',width=8,height=8)
ggsave(g_fig1, file='outputs/g_fig1.png',width=8,height=8)


results_45$cooling_effect_degc %>% min
results_45$cooling_effect_degc %>% max

results_45$delta_t_degc %>% min
results_45$delta_t_degc %>% max

lm(cooling_effect_degc~delta_t_degc, data=results_45) %>% summary
