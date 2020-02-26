list.of.packages <- c("data.table","reshape2","bsts","lubridate","dplyr","ggplot2","scales","zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only=T)

wd = "~/git/clean-cars-act-2020-analysis"
setwd(wd)

# Read in and format data
dat_march = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_March_2019.csv")
setnames(dat_march,"Count","2019-03-01")
dat_april = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_April_2019.csv")
setnames(dat_april,"Count","2019-04-01")
dat_may = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_May_2019.csv")
setnames(dat_may,"Count","2019-05-01")
dat_june = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_June_2019.csv")
setnames(dat_june,"Count","2019-06-01")
dat_august = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_August_2019.csv")
setnames(dat_august,"Count","2019-08-01")
dat_september = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_September_2019.csv")
setnames(dat_september,"Count","2019-09-01")
dat_december = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_December_2019.csv")
setnames(dat_december,"Count","2019-12-01")
dat_jan20 = fread("data/MVA_Electric_and_Hybrid_Vehicle_Registrations_by_County_as_of_January_2020.csv")
setnames(dat_jan20,"Count","2020-01-01")
dat_reg = merge(dat_march,dat_april)
dat_reg = merge(dat_reg,dat_may)
dat_reg = merge(dat_reg,dat_june)
dat_reg = merge(dat_reg,dat_august)
dat_reg = merge(dat_reg,dat_september)
dat_reg = merge(dat_reg,dat_december)
dat_reg = merge(dat_reg,dat_jan20)

dat_melt = melt(dat_reg,id.vars=c("Fuel_Category","County"),variable.name="Date")
dat_cast = dcast(dat_melt,County+Date~Fuel_Category)
names(dat_cast) = make.names(names(dat_cast))
dat_tab = data.table(dat_cast)[,.(
  PHEV=sum(Plug.In.Hybrid),
  BEV=sum(Electric)
),by=.(Date)]
dat_tab$Date = as.Date(dat_tab$Date)
missing_months = data.frame(
  Date=as.Date(
    c(
      "2019-01-01",
      "2019-02-01",
      "2019-07-01",
      "2019-10-01",
      "2019-11-01"
      )
  )
)
dat_tab = rbindlist(list(dat_tab,missing_months), fill=T)
dat_tab = dat_tab[order(dat_tab$Date)]

dat_alliance = fread("data/alliance_dat.csv")
dat_alliance$Date = as.Date(dat_alliance$Date,format="%m/%d/%y")
dat_alliance$BEV = cumsum(dat_alliance$BEV)
dat_alliance$PHEV = cumsum(dat_alliance$PHEV)

# Combine all our sources
dat = rbindlist(list(dat_alliance,dat_tab),fill=T)
dat = dat[order(dat$Date)]
dat$BEV = na.spline(dat$BEV)
dat$PHEV = na.spline(dat$PHEV)

dat$BEV = c(NA,diff(dat$BEV))
dat$PHEV = c(NA,diff(dat$PHEV))
dat = subset(dat,!is.na(BEV))

dat$ZEV = dat$BEV + dat$PHEV

zev_y = ts(dat$ZEV, frequency=12, start=c(2011,1))

### Run the bsts model
zev_ss <- AddLocalLinearTrend(list(), zev_y)
zev_ss <- AddSeasonal(zev_ss, zev_y, nseasons = 12)
zev_bsts.model <- bsts(zev_y, state.specification = zev_ss, niter = 500, ping=0, seed=2016)

### Get a suggested number of burn-ins
zev_burn <- SuggestBurn(0.1, zev_bsts.model)

# Predict until December 2024
zev_p <- predict.bsts(zev_bsts.model, horizon = 60, burn = zev_burn, quantiles = c(.025, .975))
zev_p.ts = ts(rep(NA,60),frequency=12,start=c(2020,2))
### Actual versus predicted
zev_d2 <- data.frame(
  # fitted values and predictions
  c(as.numeric(-colMeans(zev_bsts.model$one.step.prediction.errors[-(1:zev_burn),])+zev_y),  
    as.numeric(zev_p$mean)),
  # actual data and dates 
  as.numeric(c(zev_y,rep(NA,60))),
  c(as.Date(time(zev_y)),as.Date(time(zev_p.ts))))
names(zev_d2) <- c("Fitted", "Actual", "Date")

### 95% forecast credible interval
zev_posterior.interval <- cbind.data.frame(
  as.numeric(zev_p$interval[1,]),
  as.numeric(zev_p$interval[2,]), 
  subset(zev_d2, Date>as.Date("2020-01-01"))$Date)
names(zev_posterior.interval) <- c("LL", "UL", "Date")

### Join intervals to the forecast
zev_d3 <- left_join(zev_d2, zev_posterior.interval, by="Date")
zev_d3$Fitted[which(zev_d3$Fitted<0)] = 0
zev_d3$LL[which(zev_d3$LL<0)] = 0
zev_d3$UL[which(zev_d3$UL<0)] = 0

d4 = copy(zev_d3)
d4$actual_plus_estimate = d4$Actual
d4$actual_plus_LL = d4$Actual
d4$actual_plus_UL = d4$Actual
d4$actual_plus_estimate[which(is.na(d4$actual_plus_estimate))] = d4$Fitted[which(is.na(d4$actual_plus_estimate))]
d4$actual_plus_LL[which(is.na(d4$actual_plus_LL))] = d4$LL[which(is.na(d4$actual_plus_LL))]
d4$actual_plus_UL[which(is.na(d4$actual_plus_UL))] = d4$UL[which(is.na(d4$actual_plus_UL))]
d4$estimated_sum = cumsum(d4$actual_plus_estimate)
d4$LL_sum = cumsum(d4$actual_plus_LL)
d4$UL_sum = cumsum(d4$actual_plus_UL)

sum.p = ggplot(data=d4, aes(x=Date)) +
  geom_line(data=subset(d4,Date<as.Date("2020-02-01")),aes(y=estimated_sum, colour = "Actual"), size=1.2) +
  geom_line(data=subset(d4,Date>=as.Date("2020-02-01")),aes(y=estimated_sum, colour = "Estimate"), size=1.2, linetype=1) +
  theme_bw() + theme(legend.title = element_blank()) + ylab("") + xlab("") +
  geom_vline(xintercept=as.numeric(as.Date("2020-02-01")), linetype=2) + 
  geom_ribbon(data=subset(d4,Date>=as.Date("2020-02-01")),aes(ymin=LL_sum, ymax=UL_sum, fill="95% conf."), alpha=0.5) +
  scale_fill_manual("", values=c("95% conf."="grey")) +
  scale_y_continuous(labels=scales::comma) +
  expand_limits(y=300000) +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  labs(title="Cumulative Maryland ZEV sales from 2011- Jan. 2020, forecast to 2025")
ggsave("output/cumulative_projection.png", sum.p, width=10, height=5)
