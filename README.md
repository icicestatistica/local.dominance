# dominancia <img src="man/figures/logo.png" align="right" height="139" alt="" />

## Visão Geral

O pacote `local.dominance` implementa um método para detectar **regiões de dominância local** entre duas distribuições. Em vez de apenas responder "qual grupo tem média maior?", o método identifica *onde* no suporte da variável uma distribuição domina a outra.

### Principais funções

|Função|Descrição|
|-|-|
|`local_dominance()`|Função principal: roda o teste e retorna resultados + gráficos|
|`report_dominance()`|Imprime um relatório textual interpretativo|
|`plot_densities()`|Visualiza as densidades dos dois grupos|

\---

## Instalação

```r
# install.packages("devtools")
devtools::install_github("icicestatistica/local.dominance")
```

\---

## Exemplo de uso

```r
library(local.dominance)

set.seed(42)

# Grupo A: distribuição normal padrão
yA <- rnorm(200, mean = 0, sd = 1)

# Grupo B: mistura bimodal
yB <- c(rnorm(100, mean = -1, sd = 1), rnorm(100, mean = 2, sd = 1))

# Rodar o teste de dominância local
resultado <- local\_dominance(
  yA, yB,
  name_A = "Controle",
  name_B = "Tratamento",
  name_y = "Score",
  B = 999
)

# Ver resultados
resultado$p_value        # p-valor do teste omnibus
resultado$padrao         # padrão de dominância (ex: "Controle < Tratamento")
resultado$LI_assertividade  # limite inferior de assertividade (%)
resultado$partition      # tabela com as faixas dominantes

# Exibir gráfico
resultado$grafico

# Relatório completo
report_dominance(resultado)
```

\---

## Método

O método opera em três etapas:

**1. Teste omnibus (permutação)**
Estima as densidades de `yA` e `yB` via KDE e calcula `g(x) = fA(x) - fB(x)`. A estatística `D = Σ(integral\_k²)` agrega as integrais por trecho de sinal. Um teste de permutação produz o p-valor global.

**2. Identificação das regiões dominantes**
Se `p < 0.05`, compara as integrais acumuladas observadas com o quantil de referência das distribuições nulas. Changepoints são localizados por interpolação linear nos cruzamentos de `g(x) = 0`.

**3. Bootstrap de concordância**
Réplicas bootstrap estimam, em cada ponto do grid, a proporção de réplicas que concorda com o label atribuído. A assertividade final é o percentual (ponderado pela densidade conjunta) de pontos com concordância > 95%.

\---

## Interpretação da saída

* **`p_value`**: p-valor do teste omnibus. Se > 0.05, não há evidência de dominância local.
* **`padrao`**: sequência de dominância, ex.: `"A < B"` indica que B domina A em alguma região e A domina B em outra.
* **`partition`**: data frame com start/end de cada faixa, o grupo dominante, as proporções de cada grupo naquela faixa (`prop\_A`, `prop\_B`) e a razão de proporções (`RP`).
* **`LI_assertividade`**: percentual do suporte com concordância bootstrap > 95%. Quanto mais próximo de 100, mais estável é a classificação.

\---

## Licença

MIT © Isabelle Cristina Idalgo Carnielli

