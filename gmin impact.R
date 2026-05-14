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

data_slot = read_excel('nph17626-sup-0002-tabless1-s2.xlsx',sheet='Table_S1',skip=1) %>%
  rename(Tair_degC=`Temperature (C)`) %>%
  mutate(gmin_mol = Gmin_mmol / 1000)

ggplot(data_slot, aes(x=Tair_degC,y=gmin_mol)) + 
  geom_point() +
  facet_wrap(~Species) +
  scale_y_sqrt()

coef_table_slot = data_slot %>%
  group_by(Species) %>%
  #filter(Tair_degC > 30) %>%
  do(data.frame(t(coef(lm(log(gmin_mol) ~ Tair_degC, data = .))))) %>%
  rename(intercept=2,slope=3)

data_area1 = read_excel('LMA_Source_PCE2021.xlsx') %>% select(Species, Width_cm) %>% na.omit
data_area2 = read_excel('LMA extra.xlsx') %>% select(Species, Width_cm) # 3 species estimated from Kew herbarium sheets

data_slot_joined = data_slot %>% left_join(rbind(data_area1, data_area2) ,by='Species')



process_species <- function(species_this)
{
  Tair_this = data_slot_joined %>% filter(Species==species_this) %>% pull(Tair_degC)
  gmin_mol_this = data_slot_joined %>% filter(Species==species_this) %>% pull(gmin_mol)
  df_this = data.frame(Tair=Tair_this,gmin_mol=gmin_mol_this)
  
  width_cm_this = data_slot_joined %>% filter(Species==species_this) %>% pull(Width_cm) %>% mean
  
  seq_Tair = seq(30,50,by=1)
  intercept_this = coef_table_slot %>% filter(Species==species_this) %>% pull(intercept)
  slope_this = coef_table_slot %>% filter(Species==species_this) %>% pull(slope)
  seq_gs = exp(intercept_this + slope_this * seq_Tair)
  
  
  params = data.frame(Wind=1,
                      Wleaf = width_cm_this / 100, # convert from cm to m
                      StomatalRatio = 1,  # assume all hypostomatous per martijn
                      LeafAbs = 0.86,
                      gs=seq_gs,
                      Tair=seq_Tair,
                      VPD=RHtoVPD(RH=0.6,TdegC=seq_Tair,Pa=101))
  
  
  tleaf_with_gs = apply(params, 1, FUN=function(x) {
    x = as.list(x)
    do.call("FindTleaf",args=x)
    })
  tleaf_without_gs = apply(params, 1, FUN=function(x) {
    x = as.list(x)
    x$gs = 0
    do.call("FindTleaf",args=x)
    })
  
  params$tleaf_with_gs = tleaf_with_gs
  params$tleaf_without_gs = tleaf_without_gs
  params$cooling_effect_degc = params$tleaf_with_gs - params$tleaf_without_gs
  params$Tleaf_minus_Tair_degc = params$tleaf_with_gs - params$Tair
  
  # plot the gmin fit
  g1 = ggplot(df_this, aes(x=Tair,y=gmin_mol)) +
    geom_point() +
    geom_line(data=params, aes(x=Tair,y=gs),color='purple') +
    theme_bw() +
    ggtitle(species_this)
  
  g2 =ggplot(params, aes(x=Tair,y=cooling_effect_degc)) +
    geom_line() +
    theme_bw()
  
  g3 = ggplot(params, aes(x=Tair,y=Tleaf_minus_Tair_degc)) +
    geom_line() +
    theme_bw()
  
  ggsave(ggarrange(g1, g2, g3, align='hv',nrow=1,ncol=3),
         file=sprintf('outputs/result_%s.pdf',species_this),width=10,height=3.5)
  
  cat('.')
  
  return(data.frame(species=species_this,
                    delta_t = params %>% filter(Tair==45) %>% pull(Tleaf_minus_Tair_degc),
                    cooling_tdegc = params %>% filter(Tair==45) %>% pull(cooling_effect_degc)))
}


species_list = sort(unique(data_slot$Species))
# run all species
cooling_effect_45 = do.call('rbind',lapply(species_list, process_species))

g_cooling = ggplot(cooling_effect_45, aes(x=reorder(species, cooling_tdegc),y=cooling_tdegc)) +
  geom_bar(stat='identity',fill='darkblue') +
  coord_flip() +
  theme_bw() +
  ylab('CE (°C)') +
  xlab('Species') +
  geom_hline(yintercept = 0, color='black')
ggsave(g_cooling,file='outputs/g_cooling.pdf',width=7,height=7)

g_deltat = ggplot(cooling_effect_45, aes(x=reorder(species, delta_t),y=delta_t)) +
  geom_bar(stat='identity',fill='darkred') +
  coord_flip() +
  theme_bw() +
  ylab('∆T (°C)') +
  xlab('')
ggsave(g_deltat,file='outputs/g_deltat.pdf',width=7,height=7)



g_pairs = ggplot(cooling_effect_45, aes(x=delta_t, y=cooling_tdegc)) + 
  geom_point(color='black') +
  xlim(-1,10.5) + ylim(-2,0) +
  geom_hline(yintercept = 0,color='gray') +
  geom_vline(xintercept = 0,color='gray') +
  theme_bw() +
  xlab('∆T (°C)') +
  ylab('CE (°C)')
ggsave(g_pairs, file='outputs/g_pairs.pdf', width=7,height=7)


g_fig1 = ggarrange(g_cooling, g_deltat,g_pairs, labels='AUTO',nrow=2,ncol=2)
ggsave(g_fig1, file='outputs/g_fig1.pdf',width=8,height=8)
ggsave(g_fig1, file='outputs/g_fig1.png',width=8,height=8)

#\nat Tair=45°C, RH=60%, gs=0

cooling_effect_45$cooling_tdegc %>% min
cooling_effect_45$cooling_tdegc %>% max

cooling_effect_45$delta_t %>% min
cooling_effect_45$delta_t %>% max

lm(cooling_tdegc~delta_t, data=cooling_effect_45) %>% summary
