# Baseline see Hoanan Lau et al. 2019
m.hl <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ 
		RPV + Age, data = TCGA)
summary(m.hl)

# CCA projections
m.1 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		rRad1,  data=TCGA.proj)
summary(m.1)
m.2 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		rRad2,  data=TCGA.proj)
summary(m.2)
m.3 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRad3, data=TCGA.proj)
summary(m.3)
m.4 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA1, data=TCGA.proj)
summary(m.4)
m.5 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA2, data=TCGA.proj)
summary(m.5)
m.6 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA3, data=TCGA.proj)
summary(m.6)
m.cust <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age  +  Stage +
		rRNA1 + rRad2 , data=TCGA.proj)
summary(m.cust)

