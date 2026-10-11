# Load required packages
library(psych)
library(GPArotation)

# Given data
data <- data.frame(
  Math1 = c(17.18, 19, 18, 17.5, 12.3, 18.5, 15.13, 18.25, 20, 16.7),
  Math2 = c(18.75, 20, 17, 16, 15, 15, 16, 20, 19.5, 18.25),
  English1 = c(16.75, 17.5, 18, 17.5, 16.75, 17, 18, 18.75, 17, 18.75),
  English2 = c(18.5, 20, 17, 15, 18.9, 18, 16.75, 19.9, 17.5, 20),
  Programming1 = c(19, 18, 18, 14.5, 18, 14, 13, 20, 18, 20),
  Programming2 = c(18.75, 18, 17, 17, 12, 18.5, 17.25, 20, 17.5, 20)
)


# Observational Factor Analysis
# Matrix C with 2 rows and 6 columns
C <- matrix(c(0.5, 0.5, 0.5, 0, 0, 0,
              0, 0, 0, 0.5, 0.5, 0.5), byrow = TRUE, nrow = 2)

U <- scale(data)
Factors <- U %*% t(C)
factor_loads <- cor(data, Factors)
cor(Factors[, 1], Factors[, 2])

# Orthogonal Factor Analysis with psych
ofa_result <- fa(data, nfactors = 4, rotate = "varimax")
print(ofa_result)
fa.sort(ofa_result, polar = FALSE)
f <- ofa_result$scores
u <- scale(data) - f %*% t(ofa_result$loadings)

# Oblique Factor Analysis with psych
ofa_result1 <- fa(data, nfactors = 4, rotate = "oblimin")
print(ofa_result1)
fa.sort(ofa_result1, polar = FALSE)

# Maximum-Likelihood Factor Analysis
class <- factanal(data, factors = 2)
class
apply(class$loadings^2, 1, sum) # Communality
1 - apply(class$loadings^2, 1, sum) # Uniqueness
Lambda <- class$loadings
Psi <- diag(class$uniquenesses)
S <- class$correlation
Sigma <- Lambda %*% t(Lambda) + Psi
round(S - Sigma, 6)

# Interpretation of the factors with scatter plots
class_none <- factanal(data, factors = 2, rotation = "none", scores = "Bartlett")
class_varimax <- factanal(data, factors = 2, rotation = "varimax", scores = "Bartlett")
class_promax <- factanal(data, factors = 2, rotation = "promax", scores = "Bartlett")

par(mfrow = c(1, 3))
plot(class_none$loadings[, 1],
     class_none$loadings[, 2],
     xlab = "Factor 1",
     ylab = "Factor 2",
     ylim = c(-1, 1),
     xlim = c(-1, 1),
     main = "No rotation")
text(class_none$loadings[, 1] - 0.08,
     class_none$loadings[, 2] + 0.08,
     colnames(data),
     col = "blue")
abline(h = 0, v = 0)

plot(class_varimax$loadings[, 1],
     class_varimax$loadings[, 2],
     xlab = "Factor 1",
     ylab = "Factor 2",
     ylim = c(-1, 1),
     xlim = c(-1, 1),
     main = "Varimax rotation")
text(class_varimax$loadings[, 1] - 0.08,
     class_varimax$loadings[, 2] + 0.08,
     colnames(data),
     col = "blue")
abline(h = 0, v = 0)

plot(class_promax$loadings[, 1],
     class_promax$loadings[, 2],
     xlab = "Factor 1",
     ylab = "Factor 2",
     ylim = c(-1, 1),
     xlim = c(-1, 1),
     main = "Promax rotation")
text(class_promax$loadings[, 1] - 0.08,
     class_promax$loadings[, 2] + 0.08,
     colnames(data),
     col = "blue")
abline(h = 0, v = 0)
