library(maxLik)
library(np)

#####################################
### --- CS single-index model --- ###
#####################################

# --- Data simulation --- #
n <- 2000
d <- 2
x <- cbind(rnorm(n), rnorm(n))
score <-  -3 + x %*% c(1, 2)
y <- ifelse(1/(1 + exp(-score) ) > runif(n), 1, 0)
w <- rchisq(n, df = 3) * rnorm(n, mean = 7, sd = .8)*100 

# --- Initial value --- #
pre.fit <- glm(y ~ x - 1, family = binomial(link = "logit"))
st.val <- as.vector(coef(pre.fit))
start_sim <- st.val / st.val[1]

# --- Restricted optimization --- #
R <- 50
A_si <- rbind(diag(d), - diag(d))
B_si <- rep(R/d, 2*d)

# --- Model estimation --- #
cssim <- maxLik(cssim_of, 
                method = "BFGS", 
                start = start_sim,
                a_cost = .1, b_cost = 10, c_cost = 1,
                bw = length(y)**(-1/5),
                y = y,
                x = x,
                w = w, 
                lambda = 1,
                iterlim = 200,
                control=list(tol=-1, reltol=.0001, gradtol=1e-15, printLevel=0),
                constraints = list(ineqA = A_si, ineqB = B_si)
)

# --- Estimated model --- #
theta_si_lf <- c(cssim$estimate, length(y)**(-1/5))

# --- Predictions --- #
x_test <- cbind(rnorm(1000), rnorm(1000))
pred_cssim <- sim_estimator(theta_si_lf, y = y, x = x, x_test = x_test)