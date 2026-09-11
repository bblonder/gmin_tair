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
  mutate(stomatal_ratio = 1) %>% # based on personal communication from martijn, actually from his observations
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


# coef_table_combined = data_combined %>%
#   group_by(Species) %>%
#   #filter(Tair_degC > 30) %>%
#   do(data.frame(t(coef(lm(Gmin_mmol ~ Tair_degC + I(Tair_degC^2), data = .))))) %>%
#   rename(intercept=2,slope1=3, slope2=4)
# 
# # check fits
# pdf(file='outputs/g_gmin_fits.pdf')
# by(data_combined, data_combined$Species, function(x) {
#   plot(Gmin_mmol~Tair_degC,data=x,main=x$Species[1])
#   xv=20:50
#   coef_table_combined_this = coef_table_combined %>% filter(Species==x$Species[1])
#   lines(xv, coef_table_combined_this$intercept + coef_table_combined_this$slope1*xv+coef_table_combined_this$slope2*xv^2,col='red')
#   }) 
# dev.off()


process_species <- function(species_this, wind_speed_m_s, relative_humidity)
{
  # make data frame of observed values
  Tair_this = data_combined %>% filter(Species==species_this) %>% pull(Tair_degC)
  gmin_mmol_this = data_combined %>% filter(Species==species_this) %>% pull(Gmin_mmol) 
  df_this = data.frame(Tair_degC=Tair_this,gmin_mmol=gmin_mmol_this)
  
  # get species traits
  width_cm_this = data_combined %>% filter(Species==species_this) %>% pull(Width_cm) %>% mean
  stomatal_ratio_this = data_combined %>% filter(Species==species_this) %>% pull(stomatal_ratio) %>% mean
  
  # set up energy balance params
  params = data.frame(Wind=wind_speed_m_s,
                      Wleaf = width_cm_this / 100, # convert from cm to m
                      StomatalRatio = stomatal_ratio_this,
                      LeafAbs = 0.86,
                      gs=gmin_mmol_this / 1000, # convert from mmol to mol
                      Tair=Tair_this,
                      VPD=RHtoVPD(RH=relative_humidity,TdegC=Tair_this,Pa=101))
  
  
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
    geom_point() +
    theme_bw()
  
  g2 = ggplot(params, aes(x=Tair,y=delta_t_degc)) +
    geom_point() +
    theme_bw()
  
  ggsave(ggarrange(g1, g2, align='hv',nrow=1,ncol=2),
         file=sprintf('outputs/result_%s_wind_%.2f_rh_%.2f.pdf',species_this, wind_speed_m_s, relative_humidity),width=8,height=3.5)
  
  cat('.')
  
  return(data.frame(Species=species_this,
                    wind_speed_m_s = wind_speed_m_s,
                    Tair_degC = Tair_this,
                    gmin_mmol = gmin_mmol_this,
                    relative_humidity = relative_humidity,
                    width_cm = width_cm_this,
                    delta_t_degc = params %>% filter(Tair==Tair_this) %>% pull(delta_t_degc),
                    cooling_effect_degc = params %>% filter(Tair==Tair_this) %>% pull(cooling_effect_degc)))
}

# run all species
species_list_combined = sort(unique(data_combined$Species))
# main wind speed
results_base = do.call('rbind',lapply(species_list_combined, process_species, wind_speed_m_s=1, relative_humidity=0.5)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info

results_wind_low = do.call('rbind',lapply(species_list_combined, process_species, wind_speed_m_s=0.5, relative_humidity=0.5)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info
results_wind_high = do.call('rbind',lapply(species_list_combined, process_species, wind_speed_m_s=2, relative_humidity=0.5)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info

results_rh_low = do.call('rbind',lapply(species_list_combined, process_species, wind_speed_m_s=1, relative_humidity=0.25)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info
results_rh_high = do.call('rbind',lapply(species_list_combined, process_species, wind_speed_m_s=1, relative_humidity=0.75)) %>%
  left_join(data_combined %>% select(Species, Dataset) %>% unique, by='Species') # add in dataset info


g_cooling_effect = ggplot(results_base, aes(x=reorder(Species, cooling_effect_degc),y=cooling_effect_degc,shape=Dataset,color=Tair_degC)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  ylab('Cooling effect, CE (°C)') +
  xlab('Species') +
  geom_hline(yintercept = 0, color='black') +
  theme(axis.text.y=element_text(face='italic')) +
  scale_color_viridis_c(option='C',end = 0.8,limits=c(20,50),name=expression(paste('T'['air'], ' (°C)')))
ggsave(g_cooling_effect,file='outputs/g_cooling_effect.pdf',width=7,height=7)

g_delta_t = ggplot(results_base, aes(x=reorder(Species, delta_t_degc),y=delta_t_degc,shape=Dataset,color=Tair_degC)) +
  geom_point() +
  coord_flip() +
  theme_bw() +
  ylab(expression(paste('Leaf-air temperature difference, ', Delta, 'T (°C)'))) +
  xlab('') +
  geom_hline(yintercept = 0, color='black') +
  theme(axis.text.y=element_text(face='italic')) +
  scale_color_viridis_c(option='C',end = 0.8,limits=c(20,50),name=expression(paste('T'['air'], ' (°C)')))
ggsave(g_delta_t,file='outputs/g_delta_t.pdf',width=7,height=7)



results_high_temp = results_base %>% filter(Tair_degC >= 35)

g_pairs = ggplot(results_high_temp, aes(x=delta_t_degc, y=cooling_effect_degc,color=Tair_degC,shape=Dataset)) + 
  geom_point() +
  geom_hline(yintercept = 0,color='gray') +
  geom_vline(xintercept = 0,color='gray') +
  theme_bw() +
  xlab(expression(paste('Leaf-air temperature difference, ', Delta, 'T (°C)'))) +
  ylab('Cooling effect, CE (°C)') +
  stat_smooth(method='lm',se=TRUE,data=results_high_temp,aes(x=delta_t_degc, y=cooling_effect_degc),inherit.aes = FALSE,color='blue2') +
  scale_color_viridis_c(option='C',end = 0.8,limits=c(20,50),name=expression(paste('T'['air'], ' (°C)'))) +
  annotate('text', x=0.5,y=-1.75, label=expression(paste('for T'['air'] >= '35°C')),hjust=0,parse=TRUE)
ggsave(g_pairs, file='outputs/g_pairs.pdf', width=7,height=7)
ggsave(g_pairs, file='outputs/g_pairs.png', width=7,height=7)


g_fig1 = ggarrange(g_cooling_effect, g_delta_t, labels='AUTO',nrow=2,ncol=1, common.legend = TRUE,legend='bottom')
ggsave(g_fig1, file='outputs/g_fig1.pdf',width=4.5,height=8)
ggsave(g_fig1, file='outputs/g_fig1.png',width=4.5,height=8)



results_high_temp$cooling_effect_degc %>% 
  min
results_high_temp$cooling_effect_degc %>% 
  max

results_high_temp$delta_t_degc %>% 
  min
results_high_temp$delta_t_degc %>% 
  max

lm(cooling_effect_degc~delta_t_degc, data=results_high_temp) %>% summary


# check for variation due to wind speed aand humidity
results_wind = rbind(results_base, results_wind_high, results_wind_low)
results_rh = rbind(results_base, results_rh_low, results_rh_high)


# make sensitivity analyses
g_sensitivity_cooling_effect_wind = ggplot(results_wind, aes(x=Tair_degC, y=cooling_effect_degc)) +
  geom_point(alpha=0.5) + 
  theme_bw() +
  facet_wrap(~wind_speed_m_s,labeller=label_both) +
  xlab(expression(paste('T'['air'], ' (°C)'))) +
  ylab('Cooling effect, CE (°C)')


g_sensitivity_cooling_effect_rh = ggplot(results_rh, aes(x=Tair_degC, y=cooling_effect_degc)) +
  geom_point(alpha=0.5) + 
  theme_bw() +
  facet_wrap(~relative_humidity,labeller=label_both) +
  xlab(expression(paste('T'['air'], ' (°C)'))) +
  ylab('Cooling effect, CE (°C)')

ggarrange(g_sensitivity_cooling_effect_wind, g_sensitivity_cooling_effect_rh, align='hv',labels='AUTO',nrow=2,ncol=1)



g_sensitivity_delta_t_wind = ggplot(results_wind, aes(x=Tair_degC, y=delta_t_degc)) +
  geom_point(alpha=0.5) + 
  theme_bw() +
  facet_wrap(~wind_speed_m_s,labeller=label_both) +
  xlab(expression(paste('T'['air'], ' (°C)'))) +
  ylab(expression(paste('Leaf-air temperature difference, ', Delta, 'T (°C)')))


g_sensitivity_delta_t_rh = ggplot(results_rh, aes(x=Tair_degC, y=delta_t_degc)) +
  geom_point(alpha=0.5) + 
  theme_bw() +
  facet_wrap(~relative_humidity,labeller=label_both) +
  xlab(expression(paste('T'['air'], ' (°C)'))) +
  ylab(expression(paste('Leaf-air temperature difference, ', Delta, 'T (°C)')))

g_sensitivity = ggarrange(g_sensitivity_cooling_effect_wind, g_sensitivity_cooling_effect_rh, 
          g_sensitivity_delta_t_wind,g_sensitivity_delta_t_rh,
          align='hv',labels='AUTO',nrow=2,ncol=2)
ggsave(g_sensitivity, file='outputs/g_sensitivity.pdf',width=10,height=8)
ggsave(g_sensitivity, file='outputs/g_sensitivity.png',width=10,height=8)



# look at effects on leaf size
g_leafsize = ggplot(results_base, aes(x=width_cm, y=cooling_effect_degc, color=Tair_degC,shape=Dataset)) +
  geom_point() + 
  theme_bw() +
  scale_color_viridis_c(option='C',end = 0.8,limits=c(20,50),name=expression(paste('T'['air'], ' (°C)'))) +
  scale_x_sqrt() +
  ylab('Cooling effect, CE (°C)') +
  xlab('Leaf width (cm)') +
  geom_smooth(method='lm',inherit.aes = FALSE,aes(x=width_cm, y=cooling_effect_degc))
ggsave(g_leafsize, file='outputs/g_leafsize.pdf', width=6,height=6)
ggsave(g_leafsize, file='outputs/g_leafsize.png', width=6,height=6)


lm(cooling_effect_degc ~ width_cm, data=results_base) %>% summary
