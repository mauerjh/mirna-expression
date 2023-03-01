
##########
## boxplots nanosight
scale_fill_manual(labels = c("Control", "Incident", "Remittent", "Persistent"), values=RColorBrewer::brewer.pal(3,'BuPu'))

ggplot(nanosight, aes(x=group, y=concentracao_real, fill=group)) +
  geom_boxplot() +  
  scale_fill_manual(labels = c("Control", "Incident", "Remittent", "Persistent"), values=RColorBrewer::brewer.pal(4,'BuPu')) +
           ylab("Particles/mL") +
  labs(fill = "Trajectory", size = 14) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) + geom_beeswarm()

ggplot(nanosight, aes(x=group, y=tamanho_mode_average, fill=group)) +
  geom_boxplot() +  
  scale_fill_manual(labels = c("Control", "Incident", "Remittent", "Persistent"), values=RColorBrewer::brewer.pal(4,'BuPu')) +
  ylab("Mode of particle size (nm)") +
  labs(fill = "Trajectory", size = 14) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        title = element_text(size=12), legend.text=element_text(size=12)) + geom_beeswarm()
