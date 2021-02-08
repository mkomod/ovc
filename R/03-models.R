# Baseline see Hoanan Lau et al. 2019
m.hl <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		RPV, data = TCGA)
m.hl_p <- survival::coxph(Surv(Progression_free_survival_days, PFS_event) ~ Age +
		    RPV, data=TCGA.proj)


# CCA projections ---
m.1 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		rRad1,  data=TCGA.proj)
m.2 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		rRad2,  data=TCGA.proj)
m.3 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRad3, data=TCGA.proj)
m.4 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA1, data=TCGA.proj)
m.5 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA2, data=TCGA.proj)
m.6 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage + 
		rRNA3, data=TCGA.proj)
m.cust <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age  +  Stage +
		rRNA1 + rRad2 , data=TCGA.proj)


# PCA projections ---
p.1 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1,  data=TCGA.pc.proj)
p.2 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1 + pRad2,  data=TCGA.pc.proj)
p.3 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1 + pRad2 + pRad3,  data=TCGA.pc.proj)
p.4 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRNA1,  data=TCGA.pc.proj)
p.5 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRNA1 + pRNA2,  data=TCGA.pc.proj)
p.6 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRNA1 + pRNA2 + pRNA3,  data=TCGA.pc.proj)
p.7 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1 + pRNA1,  data=TCGA.pc.proj)
p.8 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1 + pRNA1 + pRad2 + pRNA2,  data=TCGA.pc.proj)
p.9 <- survival::coxph(Surv(Overall_survival_days, OS_event) ~ Age + Stage +
		pRad1 + pRNA1 + pRad2 + pRNA2 + pRad3 + pRNA3,  data=TCGA.pc.proj)

