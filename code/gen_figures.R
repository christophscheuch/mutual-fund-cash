################################################################################
# load packages -----------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)
library(tikzDevice)
library(nleqslv)

################################################################################
# set global parameters ---------------------------------------------------
omega <- 25
gamma <- 277
Fee <- 0.045
f <- 0.015
a <- 0.005
phi_0 <- 0.065

# domain for fund flow sensitivity
rspace <- seq(-0.5, 0.5, 0.001) 

# domain for default threshold
grid.lambda <- seq(from = 0.001, to = 0.1, by = 0.001)  
grid.eta <- grid.lambda

# parameters and domain for redemption shock model
lT = 0.01
p = 0.2
rho_line <- seq(0, 8, by = 0.01)

################################################################################
# graphics settings and plotting tools ------------------------------------
options(tikzDocumentDeclaration = "\\documentclass[12pt]{article}")

theme_custom <- function() {
  theme_bw() + theme(legend.title = element_blank(),
                     legend.key = element_blank(),
                     legend.background = element_blank(),
                     plot.title = element_text(hjust = 0.5))
}

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, 
                                       position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - 
                                                               lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - 
                                                             lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
}

################################################################################
# comparative statics for fund flow sensitivity ---------------------------
ffs <- function(r, phi, lambda, eta, t, default = TRUE) {
  
  sigma_q <- (omega/(gamma + (t-1)*omega)*2*a)^2/omega
  sigma_qt <- (omega/(gamma + t*omega)*2*a)^2/omega
  
  qR <- (phi-lambda/2+eta) / (2*a)
  q <- (qR*(phi+eta)-a*qR^2-lambda*sigma_q/sqrt(2*pi))/(f+eta)
  qC <- q-qR

  phi_bar <- -eta+sqrt((lambda/2)^2+4*a*(eta+f)*
                           (Fee/f+2*lambda*sigma_qt/sqrt(2*pi)/(f+eta)))
  
  # default threshold
  rstar <- ((phi+(gamma+t*omega)/omega*(phi_bar-phi))*qR-a*qR^2-
              eta*qC-f*q-lambda*(qR-(phi_bar-lambda/2+eta)/(2*a)))/q
  
  if (r < rstar & default == TRUE) {
    flow <- -1
    return(flow)
  } else {
    
    # selling treshold
    rtilde <- (phi*qR-a*qR^2-eta*qC-f*q)/q
    if (r >= rtilde) {
      R <- r*q/qR + (a*qR^2+eta*qC+f*q)/qR
      phit <- phi + omega/(gamma+t*omega)*(R-phi)
      qRt <- (phit-lambda/2+eta)/(2*a)
    } else {
      qRt <- ((gamma+(t-1)*omega)/(gamma+t*omega)*phi-lambda/2+eta+ 
                omega/(gamma+t*omega)*(q*r+(a*qR^2+eta*qC+f*q+lambda*qR))/qR)/
        (2*a+lambda/qR*omega/(gamma+t*omega))
      phit <- 2*a*qRt+lambda/2-eta
      # R <- r*q/qR + (a*qR^2+eta*qC+f*q+lambda*(qR-qRt))/qR
    }
    qt <- (qRt*(phit+eta)-a*qRt^2-lambda*sigma_qt/sqrt(2*pi))/(f+eta)
    flow <- (qt-(1+r)*q)/q
    return(flow)
  }
}
ffs <- Vectorize(ffs, "r")

# figure 2: ffs for different etas without default
eta.comp1 <- as.data.frame(cbind(rspace, 
                                 ffs(rspace, phi_0, 0, 0, 1, default = F), 
                                 ffs(rspace, phi_0, 0, 0.01, 1, default = F),
                                 ffs(rspace, phi_0, 0, 0.02, 1, default = F)))
names(eta.comp1) <- c("r", "eta1", "eta2", "eta3")

eta.comp2 <- as.data.frame(cbind(rspace, 
                                 ffs(rspace, phi_0, 0, 0, 10, default = F), 
                                 ffs(rspace, phi_0, 0, 0.01, 10, default = F),
                                 ffs(rspace, phi_0, 0, 0.02, 10, default = F)))
names(eta.comp2) <- c("r", "eta1", "eta2", "eta3")

g3 <- ggplot(eta.comp1, aes(r)) +
  geom_line(aes(y = eta1, linetype = "$\\eta = 0$"), size = 1) +
  geom_line(aes(y = eta2, linetype = "$\\eta = 0.01$"), size = 1) +
  geom_line(aes(y = eta3, linetype = "$\\eta = 0.02$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 1$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

g4 <- ggplot(eta.comp2, aes(r)) +
  geom_line(aes(y = eta1, linetype = "$\\eta = 0$"), size = 1) +
  geom_line(aes(y = eta2, linetype = "$\\eta = 0.01$"), size = 1) +
  geom_line(aes(y = eta3, linetype = "$\\eta = 0.02$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 10$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

tikz(file = "output/ffs_eta.tex", width = 6, height = 3.5)
grid_arrange_shared_legend(g3, g4, nrow = 1, ncol = 2)
dev.off()

# figure 3: ffs for different lambdas without default
lambda.comp1 <- as.data.frame(cbind(rspace, 
                                    ffs(rspace, phi_0, 0, 0, 1, default = F), 
                                    ffs(rspace, phi_0, 0.02, 0, 1, default = F),
                                    ffs(rspace, phi_0, 0.04, 0, 1, default = F)))
names(lambda.comp1) <- c("r", "lambda1", "lambda2", "lambda3")

lambda.comp2 <- as.data.frame(cbind(rspace, 
                                    ffs(rspace, phi_0, 0, 0, 10, default = F), 
                                    ffs(rspace, phi_0, 0.02, 0, 10, default = F),
                                    ffs(rspace, phi_0, 0.04, 0, 10, default = F)))
names(lambda.comp2) <- c("r", "lambda1", "lambda2", "lambda3")

g1 <- ggplot(lambda.comp1, aes(r)) +
  geom_line(aes(y = lambda1, linetype = "$\\lambda = 0$"), size = 1) +
  geom_line(aes(y = lambda2, linetype = "$\\lambda = 0.02$"), size = 1) +
  geom_line(aes(y = lambda3, linetype = "$\\lambda = 0.04$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 1$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

g2 <- ggplot(lambda.comp2, aes(r)) +
  geom_line(aes(y = lambda1, linetype = "$\\lambda = 0$"), size = 1) +
  geom_line(aes(y = lambda2, linetype = "$\\lambda = 0.02$"), size = 1) +
  geom_line(aes(y = lambda3, linetype = "$\\lambda = 0.04$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 10$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

tikz(file = "output/ffs_lambda.tex", width = 6, height = 3.5)
grid_arrange_shared_legend(g1, g2, nrow = 1, ncol = 2)
dev.off()

# figure 5: comparison with berk & green with default
bg.compare1 <- as.data.frame(cbind(rspace, 
                                   ffs(rspace, phi_0, 0.0, 0, 1, default = T), 
                                   ffs(rspace, phi_0, 0.09, 0.02, 1, default = T)))
names(bg.compare1) <- c("r", "bg", "us")

bg.compare2 <- as.data.frame(cbind(rspace, 
                                   ffs(rspace, phi_0, 0, 0, 10, default = T), 
                                   ffs(rspace, phi_0, 0.09, 0.02, 10, default = T)))
names(bg.compare2) <- c("r", "bg", "us")

g7 <- ggplot(bg.compare1, aes(r)) +
  geom_line(aes(y = bg, linetype = "$\\lambda = 0, \\eta = 0$"), size = 1) +
  geom_line(aes(y = us, linetype = "$\\lambda = 0.09, \\eta = 0.02$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 1$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

g8 <- ggplot(bg.compare2, aes(r)) +
  geom_line(aes(y = bg, linetype = "$\\lambda = 0, \\eta = 0$"), size = 1) +
  geom_line(aes(y = us, linetype = "$\\lambda = 0.09, \\eta = 0.02$"), size = 1) +
  labs(x = "$r_t$", y = "$n_t(r_t, \\varphi_{t-1})$", 
       title = "$t = 10$") +
  coord_cartesian(ylim = c(-1, 2)) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

tikz(file = "output/ffs_compare.tex", width = 6, height = 3.5)
grid_arrange_shared_legend(g7, g8, nrow = 1, ncol = 2)
dev.off()

################################################################################
# comparative statics for default threshold -------------------------------
r_star <- function(phi, lambda, eta, t) {
  sigma_q <- (omega/(gamma + t*omega)*2*a)^2/omega
  phi_bar <- -eta + sqrt((lambda/2)^2+4*a*(eta+f)*
                           (Fee/f+lambda*sigma_q/sqrt(2*pi)/(eta+f)))
  
  qR <- (phi-lambda/2+eta)/(2*a)
  qC <- (phi^2-(lambda/2)^2-2*f*(phi-lambda/2+eta))/(4*a*(eta+f))-
    (lambda*sigma_q)/(sqrt(2*pi)/(eta+f))
  q <- qR + qC
  
  rstar <- ((phi+(gamma+t*omega)/omega*(phi_bar-phi))*qR-a*qR^2-
              eta*qC-f*q-lambda*(qR-(phi_bar-lambda/2+eta)/(2*a)))/q
  return(rstar)
}
r_star <- Vectorize(r_star)

# figure 4: rstar for different lambdas vs different etas
rstar.lambda <- tibble(rstar_1 = r_star(phi_0, grid.lambda, 0, 1),
                       rstar_10 = r_star(phi_0, grid.lambda, 0, 10),
                       lambda = grid.lambda)

p1 <- ggplot(rstar.lambda, aes(x = lambda)) +
  geom_line(aes(y = rstar_1, linetype = "$t = 1$"), size = 1) +
  geom_line(aes(y = rstar_10, linetype = "$t = 10$"), size = 1) +
  labs(x = "$\\lambda$", y = "$r^*(\\varphi_{t-1})$") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

rstar.eta<- tibble(rstar_1 = r_star(phi_0, 0, grid.eta, 1),
                   rstar_10 = r_star(phi_0, 0, grid.eta, 10),
                   eta = grid.eta)

p2 <- ggplot(rstar.eta, aes(x = eta)) +
  geom_line(aes(y = rstar_1, linetype = "$t = 1$"), size = 1) +
  geom_line(aes(y = rstar_10, linetype = "$t = 10$"), size = 1) +
  labs(x = "$\\eta$", y = "$r^*(\\varphi_{t-1})$") +
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_custom()

tikz(file = "output/rstar_lambda_eta.tex", width = 6, height = 3.5)
grid_arrange_shared_legend(p1, p2, nrow = 1, ncol = 2)
dev.off()

################################################################################
# redemption shock model --------------------------------------------------
two_period <- function(rho, p, f, l = 0.09, eta = 0.01, contract = 1) {
  phi_T1 <- phi_0
  
  C <- function(q) {
    a*q^2
  }
  
  C_prime <- function(q) {
    2*a*q
  }
  
  # compensation based on total assets under management
  if (contract == 1) {
    qR_T1 <- (phi_T1 - lT + eta) / (2*a)
    qC_T1 <- (qR_T1*phi_T1 - C(qR_T1) + eta*p*rho -lT*qR_T1 - 
                f*(qR_T1-p*rho)) / (eta + f)
    
    if (rho > qC_T1) {
      fun <- function(x) {
        tmp <- numeric(2)
        tmp[1] <- (phi_T1 - (1-p)*C_prime(x[1]) - 
                     p*C_prime(x[1]- max(rho-x[2], 0)/(1-l)) - lT - f) - 
          (p/(1-l)*phi_T1 - p/(1-l)*C_prime(x[1] - max(rho-x[2],0)/(1-l)) - 
             p/(1-l)*lT - (1-p)*eta - f*(1-p+p/(1-l)))/(1-p+p/(1-l))
        tmp[2] <- -(x[1]-p*(max(rho-x[2], 0))/(1-l))*phi_T1 + 
          (1-p)*C(x[1]) + p*C(x[1] - max(rho-x[2], 0)/(1-l)) +
          (1-p)*lT*x[1] + p*lT*(x[1] - max(rho-x[2], 0)/(1-l)) + 
          eta*(1-p)*x[2] + f*(x[1] + (1-p)*x[2] - p*(max(rho-x[2], 0))/(1-l))
        tmp
      }
      qR_T1 <- nleqslv(c(10, 5), fun)$x[1]
      qC_T1 <- nleqslv(c(10, 5), fun)$x[2]
    }
  }
  
  # compensation based on illiquid holdings only
  if (contract == 2) {
    if (rho == 0) {
      qC_T1 <- rho
      qR_T1 <- (phi_T1 - lT - f) / a
    } else {
      if (rho > ((1-l)^2*eta)/(2*a*p)) {
        qC_T1 <- rho - ((1-l)^2*eta)/(2*a*p)
      } else {
        qC_T1 <- 0
      }
      fun <- function(x) {
        tmp <- numeric(1)
        tmp[1] <- -(x-p*(rho-qC_T1)/(1-l))*phi_T1 + (1-p)*C(x) + 
          p*C(x - (rho-qC_T1)/(1-l)) + (1-p)*lT*x + 
          p*lT*(x - (rho-qC_T1)/(1-l)) + eta*(1-p)*qC_T1 + 
          f*(x - p*(rho-qC_T1)/(1-l))
        tmp
      }
      qR_T1 <- nleqslv(c(10), fun)$x
    }
  }
  
  # fund default if cash holdings happen to be negative
  if (qC_T1 < 0) {
    qC_T1 <- NA
    qR_T1 <- NA
  }
  
  # calculate total fund size and manager payoff
  q_T1 <- qR_T1 + qC_T1
  if (!is.na(qR_T1) & !is.na(q_T1)) {
    if (contract == 1) {  
      if (q_T1*f < Fee) {
        qR_T1 <- NA
        qC_T1 <- NA
        q_T1 <- NA
      }
      manager_payoff <- f*q_T1
    }
    if (contract == 2) {  
      if (qR_T1*f < Fee) {
        qR_T1 <- NA
        qC_T1 <- NA
        q_T1 <- NA
      }
      manager_payoff <- f*qR_T1
    }
  } else {
    manager_payoff <- NA
  }
  
  # return results
  return(list(manager_payoff = manager_payoff,
              fund_size = q_T1,
              illiquid_holdings = qR_T1,
              cash_holdings = qC_T1))
}
two_period <- Vectorize(two_period)

# find fee for manager compensated on illiquid holdings only without any shocks
target <- two_period(0, p = p, f = f)[, 1]$manager_payoff
fee_match <- function(x) {
  target - two_period(0, p = p, f = x, contract = 2)[, 1]$manager_payoff
}
match <- nleqslv(f, fee_match)$x

# solve for different contracts
two_period_solution <- t(two_period(rho_line, p = p, f = f)) %>%
  as_tibble() %>% mutate_all(funs(as.numeric(.)))
two_period_solution <- add_column(two_period_solution, rho = rho_line)
two_period_solution2 <- t(two_period(rho_line, p = p, f = match, contract = 2)) %>%
  as_tibble() %>% mutate_all(funs(as.numeric(.)))
two_period_solution2 <- add_column(two_period_solution2, rho = rho_line)

p1 <- ggplot(NULL, aes(x = rho, y = manager_payoff)) +
  geom_line(data = two_period_solution, 
            aes(linetype = "Fee on Illiquid \\& Liquid Holdings"), size = 1) +
  geom_line(data = two_period_solution2, 
            aes(linetype = "Fee on Illiquid Holdings"), size = 1) +
  labs(x = "$\\rho$", y = "Fee", title = "Management Compensation") +
  scale_x_continuous(expand = c(0,0)) + 
  theme_custom()

p2 <- ggplot(NULL, aes(x = rho, y = fund_size)) +
  geom_line(data = two_period_solution, 
            aes(linetype = "Fee on Illiquid \\& Liquid Holdings"), size = 1) +
  geom_line(data = two_period_solution2, 
            aes(linetype = "Fee on Illiquid Only"), size = 1) +
  labs(x = "$\\rho$", y = "$q_{T-1}$", title = "Total Fund Size") +
  scale_x_continuous(expand = c(0,0)) + 
  theme_custom()

p3 <- ggplot(NULL, aes(x = rho, y = illiquid_holdings)) +
  geom_line(data = two_period_solution, 
            aes(linetype = "Fee on Illiquid \\& Liquid Holdings"), size = 1) +
  geom_line(data = two_period_solution2, 
            aes(linetype = "Fee on Illiquid Only"), size = 1) +
  labs(x = "$\\rho$", y = "$q^R_{T-1}$", title = "Illiquid Holdings") +
  scale_x_continuous(expand = c(0,0)) + 
  theme_custom()

p4 <- ggplot(NULL, aes(x = rho, y = cash_holdings)) +
  geom_line(data = two_period_solution, 
            aes(linetype = "Fee on Illiquid \\& Liquid Holdings"), size = 1) +
  geom_line(data = two_period_solution2, 
            aes(linetype = "Fee on Illiquid Only"), size = 1) +
  labs(x = "$\\rho$", y = "$q^C_{T-1}$", title = "Cash Holdings") +
  scale_x_continuous(expand = c(0,0)) + 
  theme_custom()

tikz(file = "output/redemption_shock.tex", width = 6, height = 6)
grid_arrange_shared_legend(p1, p2, p3, p4, nrow = 2, ncol = 2)
dev.off()