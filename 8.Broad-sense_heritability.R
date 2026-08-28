library(lmerTest)
data <- read.table("trait.txt", header = TRUE)
factor_vars <- c("Line", "Loc", "Rep")
data[factor_vars] <- lapply(data[factor_vars], as.factor)
calculate_heritability <- function(trait, data) {
  data <- data[!is.na(data[[trait]]), ]
  formula <- reformulate(
    c("(1|Line)", "(1|Loc)", "(1|Line:Loc)", "(1|Loc:Rep)"),
    response = trait
  )
  model <- lmer(formula, data = data, REML = TRUE)
  vc <- as.data.frame(VarCorr(model))
  sigma2 <- attr(VarCorr(model), "sc")^2  
  Vg <- vc$vcov[vc$grp == "Line"]
  Vgl <- vc$vcov[vc$grp == "Line:Loc"]
  Ve <- sigma2
  n_lines <- nlevels(data$Line)
  n_locs <- nlevels(data$Loc)
  rep_tab <- table(data$Line, data$Loc)
  n_reps_avg <- mean(rep_tab > 0)
  denominator <- Vg + Vgl / n_locs + Ve / (n_locs * n_reps_avg)
  H2 <- Vg / denominator
  se_H2 <- sqrt(
    2 * (1 - H2)^2 *
    ((Vgl / n_locs + Ve / (n_locs * n_reps_avg))^2) /
    (denominator^2 * (n_lines - 1))
  )
  data.frame(
    Trait = trait,
    nLines = n_lines,
    nLocs = n_locs,
    avgRep = n_reps_avg,
    Vg = Vg,
    Vgl = Vgl,
    Ve = Ve,
    H2 = H2,
    SE = se_H2,
    CI_low = max(0, H2 - 1.96 * se_H2),
    CI_up = min(1, H2 + 1.96 * se_H2)
  )
}
result <- calculate_heritability("Trait", data)
print(result)
