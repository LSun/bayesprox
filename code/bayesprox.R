bayesprox <- function (y, X, D, lambda, penalty = "l1", par = NULL, method = "fixpt",
                       max_iter = 2e3, fixpt.tol = 1e-5,
                       primal_init = NULL, dual_init = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(D)
  XtX <- crossprod(X)
  XtXi <- solve(XtX)
  DXtXiDt <- D %*% XtXi %*% t(D)
  DXtXiXty <- D %*% XtXi %*% t(X) %*% y
  beta0 = .lm.fit(x = X, y = y)$coef

  if (!is.null(dual_init)) {
    rho <- dual_init
  } else {rho <- rep(0, m)}

  iter_count = 0
  converged = FALSE
  obj_vec <- obj_fun(beta0, y, X, D, lambda, penalty, par)
  beta.mat <- beta0

  if (method == "fixpt") {
    c <- 1.95 / svd(DXtXiDt)$d[1]
    w <- lambda / c
    while(!converged && iter_count < max_iter) {
      beta <- beta0 - c * XtXi %*% t(D) %*% rho
      obj_val <- obj_fun(beta, y, X, D, lambda, penalty, par)
      rho.new <- fixpt.fun(rho, DXtXiDt, DXtXiXty, c, w, penalty, par)
      beta.new <- beta0 - c * XtXi %*% t(D) %*% rho.new
      obj_val.new <- obj_fun(beta.new, y, X, D, lambda, penalty, par)
      # converged <- (sum((rho.new - rho)^2) / sum(rho^2) <= fixpt.tol)
      # converged <- (sum((beta.new - beta)^2) / sum(beta^2) <= fixpt.tol)
      converged <- (log(abs(obj_val.new - obj_val)) - log(obj_val) <= log(fixpt.tol))
      # converged <- (abs(obj_val.new - obj_val) / abs(obj_val) <= 1e-10)
      rho <- rho.new
      obj_vec <- c(obj_vec, obj_val.new)
      beta.mat <- cbind.data.frame(beta.mat, beta.new)
      iter_count <- iter_count + 1
    }
  } else if (method == "iter") {
    c <- 1.95 / svd(DXtXiDt)$d[1]
    w <- lambda / c
    while(!converged && iter_count < max_iter) {
      beta <- beta0 - c * XtXi %*% t(D) %*% rho
      obj_val <- obj_fun(beta, y, X, D, lambda, penalty, par)
      gamma <- D %*% beta
      rho.new <- (rho + gamma) - prox(rho + gamma, w, penalty, par)
      beta.new <- beta0 - c * XtXi %*% t(D) %*% rho.new
      obj_val.new <- obj_fun(beta.new, y, X, D, lambda, penalty, par)
      # converged <- (sum((rho.new - rho)^2) / sum(rho^2) <= fixpt.tol)
      # converged <- (sum((beta.new - beta)^2) / sum(beta^2) <= fixpt.tol)
      converged <- ((obj_val.new - obj_val)^2 / obj_val^2 <= fixpt.tol)
      # converged <- (abs(obj_val.new - obj_val) / abs(obj_val) <= 1e-10)
      rho <- rho.new
      iter_count <- iter_count + 1
      obj_vec[iter_count] <- obj_val
      beta.mat[, iter_count] <- beta.new
    }
  }

  beta <- beta0 - c * XtXi %*% t(D) %*% rho

  obj_val <- obj_fun(beta, y, X, D, lambda, penalty, par)

  return(list(
    beta = beta,
    beta.mat = beta.mat,
    obj_val = obj_val,
    obj_vec = obj_vec,
    iter_count = iter_count,
    converged = converged))
}

fixpt.fun <- function (rho, DXtXiDt, DXtXiXty, c, w, penalty, par) {
  rho.temp <- rho - c * DXtXiDt %*% rho + DXtXiXty
  rho.new <- (rho.temp - prox(rho.temp, w, penalty, par)) * 0.95 + rho * 0.05
  return(rho.new)
}

rss_LS = function (beta, y, X) {
  return(sum((y - X %*% beta)^2))
}

obj_fun = function (beta, y, X, D, lambda, penalty, par) {
  lmloss = rss_LS(beta, y, X)/2
  reg = lambda * phi(D %*% beta, penalty, par)
  obj_val = lmloss + reg
  return(obj_val = obj_val)
}

phi = function (beta, penalty, par) {
  if (penalty == "l1") {
    return(sum(abs(beta)))
  } else {
    if (penalty == "dp") {
      return(sum(log(1 + abs(beta) / par)))
    } else {
      if (penalty == "lq") {
        return(sum(abs(beta)^par))
      } else {
        if (penalty == "l0") {
          return(sum(beta != 0))
        } else {
          stop ("invalid penalty")
        }
      }
    }
  }
}

prox = function(z, w, penalty, pen.par) {
  if (penalty == "l1") {
    prox_z = sign(z) * pmax(abs(z) - w, 0)
    return(prox_z)
  } else {
    if (penalty == "dp") {
      a = pen.par[1]
      prox_z = rep(0, length(z))
      I1 = (abs(z) > (2 * sqrt(w) - a))
      I2 = (abs(z) >= w / a)
      zI1I2 = z[I1 & I2]
      zI1I2.abs = abs(zI1I2)
      prox_zI1I2 = (zI1I2.abs - a + sqrt((zI1I2.abs + a)^2 - 4 * w)) / 2
      prox_z[I1 & I2] = prox_zI1I2 * sign(zI1I2)

      if (sqrt(w) > a) {
        zI1nI2 = z[I1 & (!I2)]
        zI1nI2.abs = abs(zI1nI2)
        prox_zI1nI2 = (zI1nI2.abs - a + sqrt((zI1nI2.abs + a)^2 - 4 * w)) / 2
        fprox_zI1nI2 = (prox_zI1nI2 - zI1nI2.abs)^2 / 2 + w * log(1 + prox_zI1nI2 / a)
        f0 = zI1nI2.abs^2 / 2
        prox_zI1nI2 = (fprox_zI1nI2 < f0) * prox_zI1nI2
        prox_z[I1 & !I2] = prox_zI1nI2 * sign(zI1nI2)
      }

      return(prox_z)

    } else {
      if (penalty == "lq") {
        q = pen.par[1]
        bwq = (2 * w * (1 - q))^(1 / (2 - q))
        hwq = bwq + w * q * bwq^(q - 1)

        # threshold rule for lq when 0 < q < 1
        prox_z = rep(0, length(z))
        I = (abs(z) > hwq)

        gamma = abs(z[I])
        get_gamma = SQUAREM::squarem(par = gamma, fixptfn = lqfixpt, z = gamma, w = w, q = q)
        gamma = get_gamma$par

        prox_z[I] = gamma * sign(z[I])

        return(prox_z)
      } else {
        if (penalty == "l0") {
          prox_z = (abs(z) > sqrt(2 * w)) * z
          return(prox_z)
        } else {
          stop ("invalid penalty")
        }
      }
    }
  }
}
