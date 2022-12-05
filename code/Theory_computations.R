
# How does Xi change w alpha? -------------------------------------------------------------

## A 3sp consumer-resource system ---------------------------------
Xi_3CR <- tibble(a=seq(0,10,0.01)) %>%
  mutate(Xi=pmap(.,get_Xi_3CR)) %>%
  unnest(Xi) %>%
  rename(alpha=a)

ggplot(Xi_3CR) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  aes(x=alpha,y=Xi) + 
  geom_line() +
  labs(x=expression(paste(alpha)), y=expression(paste(Xi))) +
  theme(legend.position = "none")

## A n sp consumer resource system ---------------------------------
Xi_nCR <- expand_grid(s_r=c(3,6,9), a_c_mean= 0.01*2^c(0:9), 
                          run=c(1:100), a_0_mean=0.2, a_0_range=0.1) %>%
  mutate(s_c=s_r/3)%>%
  mutate(a_c_range=a_c_mean/2) %>%
  mutate(A_CR=pmap(.,make_A_CR)) %>%
  mutate(Xi_0=pmap(.,get_Xi_nCR)) %>%
  unnest(Xi_0) %>%
  mutate(alpha=as_factor(a_c_mean)) 

ggplot(Xi_nCR %>% mutate(s_r = as_factor(s_r))) + 
  scale_colour_manual(values=cbPalette) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  aes(x=alpha,y=log10(Xi_0), col=s_r) + 
  geom_boxplot() + 
  labs(x=expression(paste(bar(alpha))), y=expression(paste(log[10],Xi)),
       col=expression(paste(s[r]))) +
  coord_cartesian(ylim=c(-20,0))

## A n sp food chain ---------------------------------
# Done on a cluster by running the code below 100 times
#Xi_nchain <- expand_grid(s_r=c(2, 4), s_c=c(1, 2), s_p=1, 
#a_c_mean= round(10^seq(-1, 0.5, length.out=5), 2), 
#a_p_mean= 10^seq(-1, 1, length.out=3),
#a_0_mean=0.2, a_0_range=0.1,
#run=c(1:1)) %>%
#  mutate(a_c_range=a_c_mean/2)%>% # 
#  mutate(a_p_range=a_p_mean/2)%>%#
#  mutate(A_CR = pmap(.,make_A_CR)) %>%#make A_CR
#  mutate(A_PC = pmap(.,make_A_PC)) %>%#and A_PC
#  mutate(A_chain=pmap(.,make_A_chain))%>%#put together into one A for the chain 
#  mutate(Xi_0=pmap(.,get_Xi_nchain)) %>% #get Xi
#  unnest(Xi_0) 

# How does Xi change w environmental change? -------------------------------------------------------------

## 3sp consumer-resource -----------------------------------------------
dXi_3CR <- expand_grid(r_1=seq(-0.1, 0.1, 0.05), r_2=seq(-0.1, 0.1, 0.05), 
                          r_3g=seq(-0.1, 0.1, 0.05), r_3=seq(-0.1, 0.1, 0.05),
                          alpha=c(0.3, 3.2)) %>%
  mutate(delta=(r_1+r_2)/2-r_3+r_3g) %>%
  mutate(dXi=get_dXi_3CR_de(a=alpha, r_1=r_1, r_2=r_2, 
                            r_3=r_3, r_3g=r_3g)) %>%
  mutate(alpha=as_factor(alpha))

ggplot(dXi_3CR) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  aes(x=delta, y=dXi, col=alpha) + 
  geom_line() + 
  geom_hline(yintercept = 0, lty="dotted")+
  labs(x=expression(paste(delta)),
       y=expression(paste("d", Xi[epsilon], "/d", epsilon)),
       col=expression(paste(alpha))) +
  theme(legend.position="bottom")

## 3sp. food chain -----------------------------------------------
dXi_3chain <- expand_grid(a_2 = c(0.5, 7.5), a_3 = c(0.72, 2.68), 
                          r_1=seq(-0.2, 0.2, 0.02), r_2=seq(-0.2, 0.1, 0.02), r_2g=0.1, 
                          r_3=seq(-0.2, 0.2, 0.02), r_3g=seq(-0.2, 0.1, 0.02)) %>%
  mutate(dXi=get_dXi_3Chain_de(a_2=a_2, a_3=a_3, r_1=r_1, r_2=r_2, r_3=r_3,
                               r_2g=r_2g, r_3g=r_3g)) %>%
  mutate(delta = r_1 - r_2 + r_2g) %>%
  mutate(delta2 = r_1 - r_3 + r_3g - r_2g) %>%
  mutate(a_3=as_factor(a_3))

ggplot(dXi_3chain) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill=alpha(cbPalette[1], 0.1))) +
  scale_colour_gradient2() +
  aes(x=delta, y=delta2, col=dXi, lty=a_3) + 
  geom_point() + 
  labs(x=expression(paste(delta)),
       col=expression(paste("d", Xi[epsilon], "/d", epsilon)),
       y=expression(paste(zeta))) +
  facet_grid(a_3 ~ a_2, labeller = label_bquote(rows=paste(alpha[3], "=", .(a_3)),
                                                cols=paste(alpha[2], "=", .(a_2)))) 

## n sp. consumer resource -------------------------------
dXi_nCR <- Xi_nCR %>%
  filter(a_c_mean%in%c(0.04, 2.56)) %>% #only for a selection of a's
  expand(nesting(s_r, s_c, alpha, A_CR, Xi_0), 
                                f_c_mean=c(0.75, 1, 1.25),
                                f_r_mean=c(0.75, 1, 1.25),
                                g_mean=c(0.75, 1, 1.25)) %>%
  mutate(Finv = pmap(., get_Finv_CR)) %>% #make F inverse
  unnest_wider(Finv) %>% #flatten to also get true summaries of the fs
  mutate(G    = pmap(., get_G_C)) %>%
  unnest_wider(G) %>% #flatten to also get true summaries of the gs
  mutate(A_CR = pmap(., function(Finv_cr, G_C, A_CR,...) Finv_cr%*%(G_C%*%A_CR%*%G_C))) %>% #new A, includes env. change effects, overwrites old A
  mutate(Xi   = pmap(., get_Xi_nCR)) %>% #Xi with environmental change effects
  unnest(Xi) %>% #flatten
  mutate(theta=Xi/Xi_0) %>% #compute effect on Xi
  drop_na(theta) %>% #remove all NA and Inf (not many)
  filter(!is.infinite(theta)) %>%
  mutate(delta=(f_r_mean_true-1)+(g_c_mean_true-f_c_mean_true)) #compute summary stat (delta); f=1+re, so r=(f-1)/e, so delta proportional to delta of f's

ggplot(dXi_nCR %>% mutate(s_c=as_factor(s_c))) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  aes(x=delta, y=log10(theta), col=alpha) + 
  geom_point(size=0.25, shape=1, alpha=0.2) + 
  geom_smooth(method=lm, se=F) + 
  geom_hline(yintercept = 0, lty="dotted") + 
  labs(x=expression(paste(delta)),
       y=expression(paste(log[10], "(", Xi[epsilon], "/", Xi, ")")),
       col=expression(paste(bar(alpha)))) +
  coord_cartesian(ylim=c(-5,5)) + 
  facet_grid(. ~ s_r, labeller = label_bquote(cols=paste(s[r], "=", .(s_r))))

## A n sp food chain ---------------------------------
# Done on a cluster by running the code below 100 times
#dXi_nchain <- Xi_nchain %>% 
#expand(nesting(s_r, s_c, s_p, a_c_mean, a_p_mean, A_CR, A_PC, Xi_0), 
#       f_r_mean=c(0.75, 1, 1.25), f_c_mean=c(0.75, 1, 1.25), f_p_mean=c(0.75, 1, 1.25), 
#       g_c_mean=c(0.75, 1, 1.25), g_p_mean=c(0.75, 1, 1.25)) %>%
#  mutate(G_C = pmap(., get_G_C)) %>%
#  unnest_wider(G_C) %>%
#  mutate(G_P = pmap(., get_G_P)) %>%
#  unnest_wider(G_P) %>%
#  mutate(Finv_cr = pmap(., get_Finv_CR)) %>%
#  unnest_wider(Finv_cr) %>%
#  mutate(Finv_p  = pmap(., get_Finv_P)) %>%
#  unnest_wider(Finv_p) %>%
#  mutate(Finv = pmap(., get_Finv_RCP)) %>% #collects the Fs of the resources, and consumers (in Finv_cr) and of the predators (in Finv_p)
#  mutate(A_CR = pmap(., function(G_C, A_CR,...) (G_C%*%A_CR%*%G_C))) %>% #new A_CR, includes env. change effects, overwrites old A_CR
#  mutate(A_PC = pmap(., function(G_P, A_PC,...) (A_PC%*%G_P))) %>% #new A_PC, includes env. change effects, overwrites old A_PC  
#  mutate(A_chain = pmap(.,make_A_chain))%>%#put together into one A for the chain 
#  mutate(A_chain = pmap(., function(A_chain, Finv,...) (Finv%*%A_chain))) %>% #finally multiply with F, overwrites previous A_chain
#  mutate(Xi=pmap(.,get_Xi_nchain)) %>%
#  unnest(Xi) %>%
#  mutate(theta=Xi/Xi_0) %>%
#  drop_na(theta) %>% #remove all na and inf (not many)
#  filter(!is.infinite(theta)) %>%
#  mutate(delta = (f_r_mean_true-1) + (g_c_mean_true-1) - (f_c_mean_true-1)) %>% #f=1+re, so r=(f-1)/e, so delta proportional to delta of f's
#  mutate(delta_2 = (f_r_mean_true-1) - (g_c_mean_true-1) + 
#           (g_p_mean_true-f_p_mean_true)) %>% #f=1+re, so r=(f-1)/e, so delta proportional to delta of f's
#  mutate(delta_2_factor = as_factor(round(delta_2, 1)))

#s_c_select      <- 2 #plot for a selected level of s_c
#ggplot(dXi_nchain %>% 
#         filter(s_c==s_c_select,
#                a_p_mean %in% c(1),
#                a_c_mean %in% c(0.24, 1.33), theta>0)) + 
#  theme_bw() + 
#  theme(panel.grid = element_blank(), 
#        panel.background = element_rect(fill=alpha(cbPalette[1], 0.1))) +
#  scale_colour_gradient2() +
#  aes(x=delta, y=delta_2, col=log10(theta)) + 
#  geom_point(size=0.5, alpha=0.5) + 
#  labs(x=expression(paste(delta)),
#       y=expression(paste(zeta)),
#       col=expression(paste(log[10],"(",Xi[epsilon],"/",Xi,")")),
#       shape=expression(paste(bar(alpha[2])))) +
#  facet_grid(s_r ~  a_c_mean, 
#             labeller = label_bquote(rows=paste(s[r], "=", .(s_r)),
#                                     cols=paste(bar(alpha[2]), "=", .(a_c_mean)))) +
#  theme(axis.text.x = element_text(angle = 90))

# Spatial interpretation of Xi -----------------------------------------------

## Effects of dispersal on Xi, in a 3CR -----------------------------------------------
Xi_3CR_spatial <- expand_grid(s_r=2, s_c=1, a_c_mean=c(1:40)/10, 
                    d_max=c(1e-4, 1e-3), p=100, connectivity=0.1) %>%
  mutate(n=s_r+s_c) %>% #total nr of sp
  mutate(a_c_range = a_c_mean/2) %>%
  mutate(A = pmap(., make_A_CR)) %>% #make a local A
  mutate(A_spatial = pmap(., make_block_diagonal)) %>% #make it spatial
  mutate(R_spatial = pmap(., make_R_spatial, negative_sign=c(3))) %>% #make spatial r
  mutate(D = pmap(., make_D)) %>% #add dispersal matrices
  mutate(Xi_0 = pmap(.,get_Xi_3_spatial)) %>%
  unnest(Xi_0) 

ggplot(Xi_3CR_spatial %>% mutate(dispersal = as_factor(d_max))) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  aes(x=a_c_mean, y=Xi_0, col=dispersal) + 
  geom_point() + 
  labs(x=expression(paste(alpha)), y="prop. of sites",
       col=expression(paste(d[max])))

## Effects of env. change on Xi of 3CR  ----------
dXi_3CR_spatial <- Xi_3CR_spatial %>%
  filter(a_c_mean%in%c(0.1, 3), d_max==1e-3) %>% #only for a selection of a's and d's
  expand(nesting(s_r, s_c, n, a_c_mean, A, R_spatial, D, p, Xi_0), 
         f_r_mean=c(0.75, 1, 1.25),
         f_c_mean=c(0.75, 1, 1.25),
         g_c_mean=c(0.75, 1, 1.25),
         run=c(1:5), vary=c(0, 0.1)) %>%
  mutate(Finv = pmap(., get_Finv_CR)) %>% #make F inverse
  unnest_wider(Finv) %>% #flatten to also get true summaries of the fs
  mutate(F_ = pmap(., function(Finv_cr, ...) diag(Finv_cr)^-1)) %>% #convert Finverse to F in vector format
  mutate(R_spatial= pmap(., function(F_, R_spatial,...) F_ * R_spatial)) %>% #add effects to spatial R
  mutate(G = pmap(., get_G_C)) %>% #now get effects on consumption
  unnest_wider(G) %>% #flatten to also get true summaries of the gs
  mutate(A = pmap(., function(G_C, A,...) G_C%*%A%*%G_C)) %>% #new local A, includes env. change effects, overwrites old A
  mutate(A_spatial = pmap(., make_block_diagonal)) %>% #make it spatial
  mutate(Xi = pmap(.,get_Xi_3_spatial)) %>% #get Xi
  unnest(Xi) %>%
  mutate(theta = Xi/Xi_0) %>%
  mutate(delta=(f_r_mean_true-1)+(g_c_mean_true-f_c_mean_true)) %>% #f=1+re, so r=(f-1)/e, so delta proportional to delta of f's
  drop_na(theta) %>% #remove all na and inf (not many)
  filter(!is.infinite(theta))

ggplot(dXi_3CR_spatial) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  aes(x=delta, y=log10(theta), col=as_factor(a_c_mean)) + 
  geom_point(size=0.5, shape=1, alpha=0.3) + 
  geom_smooth(method="lm", se=F) + 
  geom_hline(yintercept = 0, lty="dotted") + 
  labs(x=expression(paste(delta)),
       y=expression(paste("log(",prop[epsilon],"/prop",")")),
       col=expression(paste(alpha))) +
  coord_cartesian(ylim=c(-0.5,0.3)) + 
  facet_grid(vary ~ ., 
             labeller = label_bquote(paste(range[alpha], "=", .(vary))))

## Effects of dispersal on Xi, in a 3sp chain -----------------------------------------------
Xi_3chain_spatial <- expand_grid(s_r=1, s_c=1, s_p=1, a_c_mean=round(seq(0.5, 10, 0.5),1), 
                                 a_p_mean=c(0.72, 10), d_max=c(1e-4, 1e-3), p=100, 
                                 connectivity=0.1) %>%
  mutate(n=s_r+s_c+s_p) %>% #total nr of sp
  mutate(a_c_range = 0) %>%
  mutate(a_p_range = 0) %>%
  mutate(A_CR = pmap(., make_A_CR)) %>% #make a local A for cons-res interactions
  mutate(A_PC = pmap(., make_A_PC)) %>% #make a local A for pred-cons interactions
  mutate(A = pmap(.,make_A_chain))%>%#put together into one A for the chain 
  mutate(A_spatial = pmap(., make_block_diagonal)) %>% #make it spatial
  mutate(R_spatial = pmap(., make_R_spatial, negative_sign=c(2, 3))) %>% #make spatial r
  mutate(D = pmap(., make_D)) %>% #add dispersal matrices
  mutate(Xi_0 = pmap(.,get_Xi_3_spatial)) %>%
  unnest(Xi_0) 

ggplot(Xi_3chain_spatial) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values=cbPalette) + 
  aes(x=a_c_mean, y=Xi_0, col=as_factor(a_p_mean)) + 
  geom_point() + 
  labs(x=expression(paste(alpha[2])), y="prop. of sites", 
       col=expression(paste(alpha[3]))) + 
  facet_grid(d_max ~ ., 
             labeller = label_bquote(paste(d[max],"=",.(d_max))))

## Effects of env. change on Xi of a 3sp food chain ----------
dXi_3chain_spatial <- Xi_3chain_spatial %>%
  filter(a_c_mean%in%c(0.5, 7.5), d_max==1e-3) %>% #only for a selection of a's and d_max
  expand(nesting(s_r, s_c, s_p, n, a_c_mean, a_p_mean, A_CR, A_PC, 
                 R_spatial, D, p, Xi_0), 
         f_r_mean=c(0.75, 1, 1.25),
         f_c_mean=c(0.75, 1, 1.25),
         f_p_mean=1,
         g_c_mean=c(0.75, 1, 1.25),
         g_p_mean=1,
         run=c(1:5), vary=c(0, 0.1)) %>% #, vary=c(0, 0.1)
  mutate(Finv = pmap(., get_Finv_CR)) %>% #make F inverse for C and R
  unnest_wider(Finv) %>% #flatten to also get true summaries of the fs
  mutate(F_CR = pmap(., function(Finv_cr, ...) diag(Finv_cr)^-1)) %>% #convert Finverse to F in vector format
  mutate(Finv = pmap(.,get_Finv_P)) %>% #make F inverse for P
  unnest_wider(Finv) %>% #flatten to also get true summaries of the fs
  mutate(F_P = pmap(., function(Finv_p, ...) diag(Finv_p)^-1)) %>% #convert Finverse to F in vector format
  mutate(R_spatial= pmap(., function(F_CR, F_P, R_spatial,...) c(F_CR, F_P) * R_spatial)) %>% #add effects to spatial R
  mutate(G = pmap(., get_G_C)) %>% #now get effects on consumption by consumers
  unnest_wider(G) %>% #flatten to also get true summaries of the gs
  mutate(G = pmap(., get_G_P)) %>% #now get effects on consumption by predators
  unnest_wider(G) %>% #flatten to also get true summaries of the gs
  mutate(A_CR = pmap(., function(G_C, A_CR,...) G_C%*%A_CR%*%G_C)) %>% #new local A_CR, includes env. change effects, overwrites old A_CR
  mutate(A_PC = pmap(., function(G_P, A_PC,...) A_PC%*%G_P)) %>% #new local A_P, includes env. change effects, overwrites old A_PC. 
  mutate(A = pmap(.,make_A_chain))%>%#put together into one A for the chain 
  mutate(A_spatial = pmap(., make_block_diagonal)) %>% #make it spatial
  mutate(Xi = pmap(.,get_Xi_3_spatial)) %>% #get Xi
  unnest(Xi) %>%
  mutate(theta = Xi/Xi_0) %>%
  mutate(delta= (f_r_mean_true-1) + (g_c_mean_true-1) - (f_c_mean_true-1)) %>% #f=1+re, so r=(f-1)/e, so delta proportional to delta of f's
  mutate(delta_2= (f_r_mean_true-1) - (g_c_mean_true-1) + 
           (g_p_mean_true-f_p_mean_true)) %>% #f=1+re, so r=(f-1)/e, so delta proportional to delta of f's
  drop_na(theta) %>% #remove all na and inf (not many)
  filter(!is.infinite(theta))

ggplot(dXi_3chain_spatial%>%filter(a_c_mean==0.5)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill=alpha(cbPalette[1], 0.1))) +
  scale_colour_gradient2() +  
  aes(x=delta, y=delta_2, col=log10(theta)) + #, col=as_factor(delta_13)) + 
  geom_point(size=1, shape=16) + 
  labs(x=expression(paste(delta)),
       y=expression(paste(zeta)),
       col=expression(paste("log(",prop[epsilon],"/prop",")"))) +
  facet_grid(vary ~ a_p_mean, 
             labeller = label_bquote(cols=paste(alpha[3], "=", .(a_p_mean)),
                                     rows=paste(range[alpha], "=", .(vary)))) +
  theme(axis.text.x = element_text(angle = 90))
  