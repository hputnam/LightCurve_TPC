### Bayesian  model showing the relationship between SGD and PR by species interactions for Jamie Kerlin
### Created on 7/27/2022
### Created by Nyssa Silbiger
### https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-Source

### Load Libraries ##########
library(tidyverse)
library(tidybayes)
library(brms)
library(posterior)
library(rstan)
library(here)
library(ggthemes)
library(purrr)

##### Load data ###########
Data<-read_csv(here("output","pi_curve_extracted_rates.csv"))

## Just use a subset of the data to test the code
Data$PAR <- as.numeric(Data$Light_Value.x)
Data$Pc <- as.numeric(Data$micromol.cm2.h)
Data_sub <- Data %>% 
  filter(!Run %in% c(4,5,6)) #%>% 
  #filter(!colony_id %in% c("Mdec-D3", "Mcav-B2", "Dlab-B7", "Dlab-D4", "Dlab-F6", "Dlab-A6", "Past-A1", "Past-A3", "Past-A6"))

Data_sub_Dlab <- Data_sub %>% filter(Species %in% "Diploria labyrinthiformis")

Data_sub_Dlab_A1 <- Data_sub %>% filter(colony_id == "Dlab-A1")
#Data_sub_Past <- Data_sub %>% filter(Species == "Porites astreoides")

ggplot(Data_sub, aes(x = PAR, y=Pc))+
  geom_point()

ggplot(Data_sub_Dlab, aes(x = PAR, y=Pc))+
  geom_point()

ggplot(Data_sub_Dlab_A1, aes(PAR, y=Pc))+
  geom_point()

#set priors
prior1 <- c(set_prior("normal(0, 5)", nlpar = "Am", lb = 0),
          set_prior("normal(0, 1)", nlpar = "AQY", lb = 0),
          set_prior("normal(0, 3)", nlpar = "Rd", lb = 0))

#model 
#Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd)

fit <- brm(bf(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), Am ~ 1, AQY ~ 1, Rd ~ 1, nl = TRUE), 
              data = Data_sub_Dlab_A1, family = gaussian(),
              prior = prior1,
              chains = 4, iter = 2000, seed = 333)

#model summary
summary(fit)

#plot model fit
plot(fit)

#posterior predictive checks
pp_check(fit)

#leave-one-out cross validation (LOO) method. The LOO assesses the predictive ability of posterior distributions 
#(a little like the pp_check function). It is a good way to assess the fit of your model. 
#You should look at the elpd estimate for each model, the higher value the better the fit. 
#By adding compare = TRUE, we get a comparison already done for us at the bottom of the summary. 
#The value with an elpd of 0 should appear, that’s the model that shows the best fit to our data.
loo(fit, compare = TRUE)

model_fit <- Data_sub_Dlab %>%
    add_predicted_draws(fit) %>%  # adding the posterior distribution
    ggplot(aes(x = PAR, y = Pc)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = Data_sub_Dlab, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("Oxygen Flux") +  # latin name for red knot
    xlab("PAR µmol m-2 s-1") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85))
model_fit

#Am = Pmax
#AQY = alpha = photochemical efficiency
#Rd_Intercept = Rdark
#NEED TO ADD THESE CALCULATIONS TO THE PARAMETER LIST
#Ik=Am/AQY
#Ic=(Am*Rd)/(AQY*(sqrt(Am^2-Rd^2))))
fixef(fit)


# Extract posterior values for each parameter
samples1 <- posterior_samples(fit, "^b")
head(samples1)

# get the predicted draws from the model
pred_draws<-fit %>% 
  epred_draws(newdata = expand_grid(PAR = seq(1,800, by = 100)), 
              re_formula = NA)

# plot the fits
ggplot(pred_draws, 
       aes(x = PAR, y = .epred)) +
  stat_lineribbon() +
  geom_point(data =pred_draws, aes(x = PAR, y =.epred ) )+
  scale_fill_brewer(palette = "Reds") +
  labs(x = "PAR", y = "Rate",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

#fit many models with a for loop

names <- unique(Data_sub_Dlab$colony_id)
names 
fits <- setNames(vector("list", length(names)), names)

for (i in names) {
  fits[[i]] <- fixef(brm(bf(Pc ~ (Am*((AQY*PAR)/(sqrt(Am^2 + (AQY*PAR)^2)))-Rd), Am ~ 1, AQY ~ 1, Rd ~ 1, nl = TRUE), 
                   data = Data_sub_Dlab, family = gaussian(),
                   prior = prior1,
                   chains = 4, iter = 2000, seed = 333))
}
fits

#why does it run with the data of all samples and not specific colonies, but put the output in the correct location?



#outputs
#https://stackoverflow.com/questions/74262090/storing-model-estimates-after-running-several-models-in-a-for-loop-in-r

#samples %>%
#  spread_draws()
