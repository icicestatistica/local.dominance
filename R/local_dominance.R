#' Visualizar densidades de dois grupos
#'
#' Gera um gráfico de densidades sobrepostas para dois grupos, útil para
#' inspecionar visualmente as distribuições antes de aplicar o teste de
#' dominância local.
#'
#' @param yA Vetor numérico com os dados do grupo A.
#' @param yB Vetor numérico com os dados do grupo B.
#' @param name_A Rótulo do grupo A. Padrão: `"A"`.
#' @param name_B Rótulo do grupo B. Padrão: `"B"`.
#' @param name_y Rótulo do eixo x. Padrão: `"y"`.
#'
#' @return Um objeto `ggplot`.
#'
#' @examples
#' set.seed(1)
#' plot_densities(rnorm(100), rnorm(100, 1), name_A = "Controle", name_B = "Tratamento")
#'
#' @export
plot_densities <- function(yA, yB, name_A = "A", name_B = "B", name_y = "y") {
  data.frame(
    y = c(yA, yB),
    x = rep(c("A", "B"), c(length(yA), length(yB)))
  ) |>
    ggplot2::ggplot(ggplot2::aes(x = y, fill = x)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::theme_minimal() +
    ggplot2::labs(fill = NULL, x = name_y, y = "Densidade") +
    ggplot2::scale_fill_manual(
      values = c("A" = "#F8766D", "B" = "#00BFC4"),
      labels = c("A" = name_A, "B" = name_B)
    )
}

# -------------------------------------------------------------------------
# Parte I — Funções internas do teste omnibus
# -------------------------------------------------------------------------

#' @keywords internal
.fill_zero_signs <- function(g) {
  s <- sign(g)
  
  zero_idx <- which(s == 0)
  if (length(zero_idx) == 0) return(g)
  
  nonzero_idx <- which(s != 0)
  if (length(nonzero_idx) == 0) return(g)
  
  for (i in zero_idx) {
    left  <- nonzero_idx[nonzero_idx < i]
    right <- nonzero_idx[nonzero_idx > i]
    
    has_left  <- length(left)  > 0
    has_right <- length(right) > 0
    
    if (has_left && has_right) {
      mid  <- (left[length(left)] + right[1]) / 2
      g[i] <- if (i <= mid) g[left[length(left)]] else g[right[1]]
    } else if (has_left) {
      g[i] <- g[left[length(left)]]
    } else {
      g[i] <- g[right[1]]
    }
  }
  g
}

#' @keywords internal
.run_integrals <- function(g, dx) {
  signs  <- sign(g)
  run_id <- cumsum(c(1L, diff(signs) != 0L))
  list(
    integrals = tapply(g * dx, run_id, sum),
    n_grid    = tabulate(run_id),
    run_id    = run_id
  )
}

#' @keywords internal
.D_stat <- function(ri) sum(ri^2)

#' @keywords internal
.any_local_dominance <- function(a, b, n_perm = 999) {
  
  pool   <- c(a, b)
  x_grid <- stats::density(pool)$x
  n_grid <- length(x_grid)
  na     <- length(a)
  nb     <- length(b)
  
  kde <- function(x) {
    d <- stats::density(x, n = n_grid)
    stats::approx(d$x, d$y, xout = x_grid, yleft = 0, yright = 0)$y
  }
  
  dx           <- x_grid[2] - x_grid[1]
  fA           <- kde(a)
  fB           <- kde(b)
  g_obs_filled <- .fill_zero_signs(fA - fB)
  
  ri_obs <- .run_integrals(g_obs_filled, dx)
  D_obs  <- .D_stat(ri_obs$integrals)
  
  ri_null <- vector("list", n_perm)
  D_null  <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    perm       <- sample(pool)
    g_p_filled <- .fill_zero_signs(kde(perm[seq_len(na)]) - kde(perm[seq_len(nb) + na]))
    ri         <- .run_integrals(g_p_filled, dx)
    ri_null[[i]] <- ri
    D_null[i]    <- .D_stat(ri$integrals)
  }
  
  list(
    x_grid  = x_grid,
    g_obs   = g_obs_filled,
    p_value = mean(D_null >= D_obs),
    D_obs   = D_obs,
    D_null  = D_null,
    fA      = fA,
    fB      = fB,
    ri_obs  = ri_obs,
    ri_null = ri_null
  )
}

# -------------------------------------------------------------------------
# Parte II — Identificação de regiões dominantes
# -------------------------------------------------------------------------

#' @keywords internal
.find_dominance_candidates <- function(out, percent_exp = 0.95) {
  
  B       <- length(out$ri_null)
  p_value <- out$p_value
  ri_obs  <- out$ri_obs
  
  order_and_acumulate <- function(x) {
    ordem     <- order(abs(x$integrals))
    integrals <- ((x$integrals)^2)[ordem]
    res    <- unlist(rep(cumsum(integrals), x$n_grid[ordem]))
    id_int <- unlist(rep(names(integrals),  x$n_grid[ordem]))
    list(acum = res, id = id_int)
  }
  
  acumulada_observado <- order_and_acumulate(ri_obs)
  observ <- acumulada_observado$acum
  id     <- acumulada_observado$id
  
  matrix_acumuladas_null <- do.call(rbind, lapply(seq_len(B), function(i) {
    order_and_acumulate(out$ri_null[[i]])$acum
  }))
  
  estimate_acumuladas <- apply(matrix_acumuladas_null, 2,
                               function(x) stats::quantile(x, percent_exp))
  
  get_dominance <- data.frame(estimate_acumuladas, observ, id) |>
    dplyr::group_by(id) |>
    dplyr::mutate(dom = dplyr::cumany(estimate_acumuladas < observ)) |>
    dplyr::group_by(id, dom) |>
    dplyr::count() |>
    dplyr::filter(dom)
  
  dominantes <- get_dominance$id
  
  dom_A_runs <- as.integer(dominantes[ri_obs$integrals[dominantes] > 0])
  dom_B_runs <- as.integer(dominantes[ri_obs$integrals[dominantes] < 0])
  
  complementar_forcado <- FALSE
  
  if (p_value > 0.05) {
    
    label <- rep("none", length(out$x_grid))
    
  } else {
    
    label <- dplyr::case_when(
      ri_obs$run_id %in% dom_A_runs ~ "A",
      ri_obs$run_id %in% dom_B_runs ~ "B",
      TRUE                          ~ "none"
    )
    
    tem_A <- length(dom_A_runs) > 0
    tem_B <- length(dom_B_runs) > 0
    
    if (!(tem_A && tem_B)) {
      
      if (p_value >= 0.01) {
        
        label <- rep("none", length(out$x_grid))
        
      } else {
        
        fila_ids   <- order(abs(ri_obs$integrals))
        fila_nomes <- names(ri_obs$integrals)[fila_ids]
        sinal_fila <- ifelse(ri_obs$integrals[fila_nomes] > 0, "A", "B")
        lado_falta <- if (!tem_A) "A" else "B"
        
        pos_complementar <- which(sinal_fila == lado_falta)[1]
        
        if (!is.na(pos_complementar)) {
          runs_adicionar <- as.integer(fila_nomes[pos_complementar:length(fila_nomes)])
          
          dom_A_runs <- unique(c(
            dom_A_runs,
            runs_adicionar[ri_obs$integrals[as.character(runs_adicionar)] > 0]
          ))
          dom_B_runs <- unique(c(
            dom_B_runs,
            runs_adicionar[ri_obs$integrals[as.character(runs_adicionar)] < 0]
          ))
          
          label <- dplyr::case_when(
            ri_obs$run_id %in% dom_A_runs ~ "A",
            ri_obs$run_id %in% dom_B_runs ~ "B",
            TRUE                          ~ "none"
          )
          
          complementar_forcado <- TRUE
        }
      }
    }
  }
  
  list(
    df                   = data.frame(x = out$x_grid, g = out$g_obs, label = label),
    complementar_forcado = complementar_forcado
  )
}

#' @keywords internal
.find_transition_regions <- function(df) {
  
  label    <- df$label
  x        <- df$x
  block_id <- cumsum(c(1L, label[-1] != label[-length(label)]))
  
  blocks <- data.frame(
    block     = unique(block_id),
    label     = tapply(label, block_id, `[`, 1),
    start_idx = tapply(seq_along(label), block_id, min),
    end_idx   = tapply(seq_along(label), block_id, max)
  )
  
  sig_blocks <- blocks[blocks$label != "none", ]
  if (nrow(sig_blocks) < 2) return(NULL)
  
  transition_regions <- list()
  
  for (i in seq_len(nrow(sig_blocks) - 1)) {
    b1 <- sig_blocks[i, ]
    b2 <- sig_blocks[i + 1, ]
    if (b1$label == b2$label) next
    
    transition_regions[[length(transition_regions) + 1]] <- data.frame(
      from      = b1$label,
      to        = b2$label,
      start     = x[b1$end_idx],
      end       = x[b2$start_idx],
      start_idx = b1$end_idx,
      end_idx   = b2$start_idx
    )
  }
  
  if (length(transition_regions) == 0) return(NULL)
  do.call(rbind, transition_regions)
}

#' @keywords internal
.find_changepoints <- function(df_final, transition_regions) {
  
  if (is.null(transition_regions) || nrow(transition_regions) == 0) return(NULL)
  
  x <- df_final$x
  g <- df_final$g
  
  cps <- lapply(seq_len(nrow(transition_regions)), function(i) {
    tr    <- transition_regions[i, ]
    idx   <- tr$start_idx:tr$end_idx
    x_sub <- x[idx]
    g_sub <- g[idx]
    
    crossings <- which(diff(sign(g_sub)) != 0)
    
    if (length(crossings) == 0) {
      x_cp <- x_sub[which.min(abs(g_sub))]
    } else {
      best <- crossings[which.min(abs(g_sub[crossings]) + abs(g_sub[crossings + 1]))]
      x_cp <- stats::approx(
        x    = g_sub[c(best, best + 1)],
        y    = x_sub[c(best, best + 1)],
        xout = 0
      )$y
    }
    
    data.frame(from = tr$from, to = tr$to, changepoint = x_cp)
  })
  
  do.call(rbind, cps)
}

#' @keywords internal
.build_partition <- function(grid_x, crossings, yA, yB) {
  
  if (is.null(crossings) || nrow(crossings) == 0) return(NULL)
  
  crossings <- crossings[order(crossings$changepoint), ]
  breaks    <- c(-Inf, crossings$changepoint, Inf)
  n_faixas  <- length(breaks) - 1
  doms      <- c(crossings$from[1], crossings$to)
  
  faixas <- lapply(seq_len(n_faixas), function(i) {
    lo  <- breaks[i]
    hi  <- breaks[i + 1]
    idx <- which(grid_x > lo & grid_x <= hi)
    if (length(idx) == 0) return(NULL)
    
    prop_A <- mean(yA > lo & yA < hi)
    prop_B <- mean(yB > lo & yB < hi)
    
    data.frame(
      faixa     = i,
      start     = ifelse(is.infinite(lo), min(grid_x), lo),
      end       = ifelse(is.infinite(hi), max(grid_x), hi),
      dominante = doms[i],
      prop_A    = round(prop_A, 4),
      prop_B    = round(prop_B, 4),
      RP        = round(prop_A / prop_B, 4)
    )
  })
  
  do.call(rbind, Filter(Negate(is.null), faixas))
}

#' @keywords internal
.estimate_labels <- function(out, partition) {
  
  df <- data.frame(x = out$x_grid, g = out$g_obs, label = "none")
  
  if (!is.null(partition)) {
    cuts      <- c(-Inf, partition$end)
    df$label  <- as.character(cut(df$x, breaks = cuts, labels = partition$dominante))
  }
  
  df
}

#' @keywords internal
.grafico_dominancia <- function(df_final, crossings, name_A = "A", name_B = "B", name_y = "y") {
  ggplot2::ggplot() +
    ggplot2::geom_point(data = df_final, ggplot2::aes(x, g, color = label)) +
    ggplot2::geom_vline(xintercept = crossings$changepoint, linetype = "dashed", color = "gray30") +
    ggplot2::geom_hline(yintercept = 0, color = "gray60", linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(color = "Dominância", x = name_y) +
    ggplot2::scale_color_manual(
      values = c("A" = "#F8766D", "B" = "#00BFC4", "none" = "gray"),
      labels = c("A" = name_A, "B" = name_B, "none" = "Sem Dominância")
    )
}

# -------------------------------------------------------------------------
# Parte III — Bootstrap de concordância
# -------------------------------------------------------------------------

#' @keywords internal
.bootstrap_concordance <- function(yA, yB, df_final, ald, B = 999) {
  
  n_grid <- nrow(df_final)
  x_grid <- df_final$x
  
  kde <- function(x, bw, xout) {
    d <- stats::density(x, bw = bw)
    stats::approx(d$x, d$y, xout = xout, yleft = 0, yright = 0)$y
  }
  
  matrix_g <- matrix(NA_real_, nrow = B, ncol = n_grid)
  
  for (b in seq_len(B)) {
    yA_b         <- sample(yA, replace = TRUE)
    yB_b         <- sample(yB, replace = TRUE)
    matrix_g[b, ] <- .fill_zero_signs(
      kde(yA_b, stats::bw.nrd0(yA_b), x_grid) -
        kde(yB_b, stats::bw.nrd0(yB_b), x_grid)
    )
  }
  
  conc_A <- colMeans(matrix_g >= 0)
  conc_B <- colMeans(matrix_g <= 0)
  
  if (all(df_final$label == "none")) {
    conc          <- pmax(conc_A, conc_B)
    expected_sign <- rep("empate", n_grid)
  } else {
    conc <- rep(NA_real_, n_grid)
    conc[df_final$label == "A"]    <- conc_A[df_final$label == "A"]
    conc[df_final$label == "B"]    <- conc_B[df_final$label == "B"]
    conc[df_final$label == "none"] <- pmax(
      conc_A[df_final$label == "none"],
      conc_B[df_final$label == "none"]
    )
    expected_sign <- df_final$label
  }
  
  peso  <- ald$fA + ald$fB
  peso  <- peso / sum(peso)
  
  acima_limiar      <- conc > 0.95
  certeza_grid      <- mean(acima_limiar, na.rm = TRUE)
  certeza_ponderada <- sum(peso * acima_limiar, na.rm = TRUE)
  
  list(
    result            = data.frame(x = x_grid, concordancia = conc, peso = peso, expected = expected_sign),
    certeza           = certeza_grid,
    certeza_grid      = certeza_grid,
    certeza_ponderada = certeza_ponderada
  )
}

#' @keywords internal
.plot_confianca <- function(bc, name_A = "A", name_B = "B", name_y = "y") {
  ggplot2::ggplot() +
    ggplot2::geom_line(data = bc$result, ggplot2::aes(x, concordancia)) +
    ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    ggplot2::labs(
      y        = "Assertividade",
      x        = name_y,
      subtitle = paste0(
        "Limite inferior para assertividade: Bruta=", round(bc$certeza, 3),
        "; Ponderada=", round(bc$certeza_ponderada, 3)
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_manual(
      values = c("A" = "#F8766D", "B" = "#00BFC4", "none" = "gray"),
      labels = c("A" = name_A, "B" = name_B, "none" = "Sem Dominância")
    )
}

# -------------------------------------------------------------------------
# Wrapper interno
# -------------------------------------------------------------------------

#' @keywords internal
.inner_local_dominance <- function(yA, yB, name_A, name_B, name_y, B, percent_exp, graficos) {
  
  ald     <- .any_local_dominance(yA, yB, n_perm = B)
  limites <- c(min(ald$x_grid), max(ald$x_grid))
  
  complementar_forcado <- FALSE
  
  if (ald$p_value < 0.05) {
    
    dc_result            <- .find_dominance_candidates(ald, percent_exp = percent_exp)
    df                   <- dc_result$df
    complementar_forcado <- dc_result$complementar_forcado
    
    if (all(df$label == "none")) {
      crossings <- data.frame(from = character(), to = character(), changepoint = numeric())
      partition <- NULL
    } else {
      trans_reg <- .find_transition_regions(df)
      crossings <- .find_changepoints(df, trans_reg)
      partition <- .build_partition(ald$x_grid, crossings, yA, yB)
    }
    
  } else {
    crossings <- data.frame(from = character(), to = character(), changepoint = numeric())
    partition <- NULL
  }
  
  df_final <- .estimate_labels(ald, partition)
  bc       <- .bootstrap_concordance(yA, yB, df_final, ald, B)
  
  if (graficos) {
    g1   <- plot_densities(yA, yB, name_A, name_B, name_y) +
      ggplot2::scale_x_continuous(limits = limites) +
      ggplot2::labs(subtitle = "Estimação das densidades")
    g2   <- .grafico_dominancia(df_final, crossings, name_A, name_B, name_y) +
      ggplot2::scale_x_continuous(limits = limites) +
      ggplot2::labs(subtitle = paste0("p-valor para o teste omnibus=", scales::pvalue(ald$p_value)))
    g3   <- .plot_confianca(bc, name_A, name_B, name_y) +
      ggplot2::scale_x_continuous(limits = limites)
    graf <- patchwork::wrap_plots(list(g1, g2, g3), ncol = 1,
                                  guides = "collect", axis_titles = "collect")
  } else {
    graf <- NULL
  }
  
  list(
    ald                  = ald,
    crossings            = crossings,
    partition            = partition,
    df_final             = df_final,
    bc                   = bc,
    graf                 = graf,
    complementar_forcado = complementar_forcado
  )
}

# -------------------------------------------------------------------------
# Funções exportadas
# -------------------------------------------------------------------------

#' Teste de Dominância Local entre Duas Distribuições
#'
#' Compara duas amostras usando estimação de densidade kernel (KDE) e um
#' teste de permutação omnibus. Identifica regiões do suporte onde uma
#' distribuição domina a outra e estima a assertividade dessa classificação
#' via bootstrap.
#'
#' @details
#' O método opera em três etapas:
#'
#' **Parte I — Teste omnibus:** estima as densidades de `yA` e `yB` via KDE
#' e calcula a diferença `g(x) = fA(x) - fB(x)` sobre um grid. A estatística
#' `D = sum(integral_k^2)` agrega as integrais por trecho de sinal. Um teste
#' de permutação produz o p-valor.
#'
#' **Parte II — Identificação de regiões dominantes:** se `p < 0.05`,
#' compara as integrais acumuladas observadas com o quantil `percent_exp` das
#' distribuições nulas para selecionar os trechos candidatos. Changepoints são
#' encontrados por interpolação linear nos cruzamentos de `g(x) = 0`.
#'
#' **Parte III — Bootstrap de concordância:** reamostras bootstrap estimam, em
#' cada ponto do grid, a proporção de réplicas que concorda com o label
#' atribuído. A assertividade final é o percentual ponderado de pontos com
#' concordância > 95%.
#'
#' @param yA Vetor numérico com os dados do grupo A.
#' @param yB Vetor numérico com os dados do grupo B.
#' @param name_A Rótulo do grupo A nos gráficos e na saída. Padrão: `"A"`.
#' @param name_B Rótulo do grupo B nos gráficos e na saída. Padrão: `"B"`.
#' @param name_y Rótulo do eixo da variável de interesse. Padrão: `"y"`.
#' @param B Número de permutações (teste omnibus) e de réplicas bootstrap
#'   (concordância). Padrão: `999`.
#' @param percent_exp Percentil das distribuições nulas usado como limiar para
#'   selecionar candidatos dominantes. Padrão: `0.95`.
#' @param graficos Lógico. Se `TRUE` (padrão), retorna um objeto `patchwork`
#'   com três painéis: densidades, dominância e assertividade.
#'
#' @return Lista com os seguintes elementos:
#' \describe{
#'   \item{estat_D}{Valor observado da estatística D.}
#'   \item{p_value}{P-valor do teste de permutação omnibus.}
#'   \item{grafico}{Objeto `patchwork` (ou `NULL` se `graficos = FALSE`).}
#'   \item{partition}{Data frame com as faixas de dominância: `start`, `end`,
#'     `dominante`, `prop_A`, `prop_B`, `RP` (razão de proporções).}
#'   \item{padrao}{String descrevendo o padrão de dominância (ex.: `"A < B"`).}
#'   \item{LI_assertividade}{Limite inferior de assertividade em percentual.}
#' }
#'
#' @examples
#' set.seed(42)
#' yA <- rnorm(200, mean = 0, sd = 1)
#' yB <- c(rnorm(100, mean = -1, sd = 1), rnorm(100, mean = 2, sd = 1))
#'
#' resultado <- local_dominance(yA, yB, name_A = "Controle", name_B = "Tratamento")
#' resultado$p_value
#' resultado$padrao
#' resultado$grafico
#'
#' @seealso [report_dominance()] para imprimir um relatório formatado.
#'
#' @export
local_dominance <- function(yA, yB, name_A = "A", name_B = "B", name_y = "y",
                            B = 999, percent_exp = 0.95, graficos = TRUE) {
  
  ld   <- .inner_local_dominance(yA, yB, name_A, name_B, name_y, B, percent_exp, graficos)
  part <- ld$partition
  part$dominante <- c("A" = name_A, "B" = name_B)[part$dominante]
  
  list(
    estat_D          = ld$ald$D_obs,
    p_value          = ld$ald$p_value,
    grafico          = ld$graf,
    partition        = part,
    padrao           = paste0(part$dominante, collapse = " < "),
    LI_assertividade = round(100 * ld$bc$certeza, 2)
  )
}

#' Relatório textual de dominância local
#'
#' Imprime no console um resumo interpretativo do resultado de
#' [local_dominance()], incluindo a estatística D, o p-valor, o padrão de
#' dominância encontrado e os gráficos gerados.
#'
#' @param result Lista retornada por [local_dominance()].
#' @param digits_D Número de casas decimais para a estatística D. Padrão: `3`.
#' @param digits_p Número de casas decimais para o p-valor. Padrão: `3`.
#' @param digits_prop Número de casas decimais para proporções. Padrão: `2`.
#'
#' @return Invisível (`NULL`). Efeito colateral: imprime texto e gráficos.
#'
#' @examples
#' set.seed(42)
#' yA <- rnorm(200, mean = 0, sd = 1)
#' yB <- c(rnorm(100, mean = -1, sd = 1), rnorm(100, mean = 2, sd = 1))
#' resultado <- local_dominance(yA, yB)
#' report_dominance(resultado)
#'
#' @export
report_dominance <- function(result, digits_D = 3, digits_p = 3, digits_prop = 2) {
  
  estat_D    <- result$estat_D
  p_value    <- result$p_value
  grafico    <- result$grafico
  partition  <- result$partition
  padrao     <- result$padrao
  LI_certeza <- result$LI_assertividade
  
  p_txt <- ifelse(
    p_value < 10^-digits_p,
    paste0("< ", format(10^-digits_p, scientific = FALSE)),
    round(p_value, digits_p)
  )
  
  if (p_value > 0.05) {
    
    cat(paste0(
      "Não encontramos transição de dominância ",
      "(t_dom = ", round(estat_D, digits_D),
      ", p-valor = ", p_txt, ").\n\n",
      "Veja o gráfico a seguir:\n\n"
    ))
    print(grafico)
    return(invisible(NULL))
    
  }
  
  cat(paste0(
    "Encontramos transição de dominância com ",
    "t_dom = ", round(estat_D, digits_D),
    " e p-valor = ", p_txt, ".\n\n",
    "O padrão encontrado foi de ", padrao,
    ", com mais de 95% de concordância em ",
    round(as.numeric(LI_certeza), digits_prop),
    "% dos pontos do domínio.\n\n",
    "As proporções encontradas estão na tabela abaixo:\n\n"
  ))
  
  print(knitr::kable(partition, digits = 3))
  cat("\n\nComo podemos ver no gráfico a seguir:\n\n")
  print(grafico)
}