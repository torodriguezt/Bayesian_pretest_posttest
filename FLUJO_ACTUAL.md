# FLUJO COMPLETO DEL PROYECTO (Paso a Paso)

## 1. EL MOTOR: BivBetaBinom.cpp

**Qué es:** Código C++ compilado con Rcpp. Contiene todas las funciones matemáticas aceleradas.

**Funciones principales que exporta:**

- `bb_constants(n1, n2, x1, x2, a0, a1, a2)` 
  - Entrada: tamaños muestrales (n1, n2), conteos (x1, x2), hiperparámetros (a0, a1, a2)
  - Salida: constantes precomputadas de la posterior (log-C, exponentes)
  - Velocidad: ~1ms
  - Costo: evita recalcular lo mismo cientos de veces

- `densBB_cpp(theta1, theta2, consts)`
  - Entrada: parámetros (θ1, θ2) y constantes precomputadas
  - Salida: valor de densidad posterior en ese punto
  - Usa: forma cerrada de la posterior (NO MCMC)

- `find_sup_H(consts)` 
  - Entrada: constantes precomputadas
  - Salida: supremo de la densidad bajo H (donde θ1=θ2)
  - Método: barrido de grilla 1D + refina
  - Necesario para: calcular e-valor del FBST

- `ev_quad(consts, sup_H, ngrid=401)`
  - Entrada: constantes, supremo bajo H, resolución de grilla
  - Salida: e-valor (un número entre 0 y 1)
  - Método: cuadratura Simpson 2D sobre todo el dominio
  - **ES LO IMPORTANTE:** Calcula integral doble SIN MCMC
  - Velocidad: ~50ms para grilla 401×401

- `simulate_evs_H(n1, n2, a0, a1, a2, M=1000)`
  - Entrada: tamaños, hiperparámetros, número de simulaciones M
  - Salida: vector de M e-valores
  - Proceso:
    1. Muestrea M valores de θ desde la PRIOR bajo H
    2. Para cada θ, simula datos binomiales (X1, X2)
    3. Calcula e-valor de esos datos
    4. Devuelve los M e-valores
  - Uso: Construir curva de error α(k) = P(ev ≤ k | H es verdadera)

- `simulate_evs_A(n1, n2, a0, a1, a2, M=1000)`
  - Igual que anterior pero desde la PRIOR completa (θ1, θ2 sin restricción)
  - Uso: Construir curva de error β(k) = P(ev > k | A es verdadera)

- `simulate_evs_H_post(n1, n2, x1_obs, x2_obs, a0, a1, a2, M=1000)`
  - Igual que `simulate_evs_H` pero muestrea desde la POSTERIOR (dados los datos observados)
  - Más realista: refleja lo que pasaría si repitieras el experimento con ese dataset

- `find_kstar(ev_H, ev_A, a=1, b=1)`
  - Entrada: vectores de e-valores bajo H y bajo A
  - Salida: k* = argmin(a·α(k) + b·β(k))
  - Significa: el cutoff que minimiza el error total
  - Rápido: O(M log M) porque enumera todos los e-valores observados

---

## 2. PIPELINE A: GENERAR TABLA A (k* genérica, prior-based)

**Archivo:** `build_article_tables.R` (líneas 76-100)

**Objetivo:** Crear una tabla de referencia: para cualquier futuro experimento con tamaños (n1, n2), cuál debería ser el cutoff k*.

**Pasos:**

```
1. sourceCpp("BivBetaBinom.cpp")
   → Compila y carga todas las funciones C++ en memoria de R

2. Define hiperparámetros KL-óptimos:
   a0 = 0.8373879, a1 = 0.8410984, a2 = 0.8053298
   → Estos minimizan divergencia Kullback-Leibler

3. Para cada combinación (n1, n2) en grilla 9×9:
   n_grid = {10, 20, 30, 40, 50, 75, 100, 150, 200}
   
   Entonces: 9×9 = 81 combinaciones
   
   Para cada una:
   ├─ ev_H ← simulate_evs_H(n1, n2, M=2000)
   │  └─ Simula 2000 datasets bajo H, calcula e-valores
   │
   ├─ ev_A ← simulate_evs_A(n1, n2, M=2000)
   │  └─ Simula 2000 datasets bajo A, calcula e-valores
   │
   └─ k_star ← find_kstar(ev_H, ev_A)
      └─ Encuentra el k que minimiza α + β
   
   → Guarda k* en matriz 9×9

4. Tiempo total: ~3 horas (por eso M=2000, no 500)

5. Salida: output/kstar_prior_table.tex
   → Tabla lista para \input{} en article.tex
```

**Por qué es importante:** Es una tabla genérica, NO depende de los datos. Cualquier investigador la puede usar para futuros estudios.

---

## 3. PIPELINE B: GENERAR TABLA B (k* por tratamiento, posterior-based)

**Archivo:** `build_article_tables.R` (líneas 102-178)

**Objetivo:** Para los 4 grupos específicos del THKS, calcular k* usando la posterior real (no la prior).

**Pasos:**

```
1. sourceCpp("BivBetaBinom.cpp")

2. Define hiperparámetros KL-óptimos (igual que Tabla A)

3. Carga datos reales del TVSFP (paquete ALA):
   datos1 <- tvsfp
   
4. Prepara 4 grupos:
   ├─ yy = "CC + TV"       (school 404)
   ├─ yn = "CC, no TV"     (school 408)
   ├─ ny = "no CC, TV"     (global)
   └─ nn = "no CC, no TV"  (school 409)

5. Para CADA grupo:
   ├─ Extrae n1, n2, x1, x2 del dataset
   │  (x1 = número de students con THKS ≥ 3 en pretest)
   │  (x2 = número de students con THKS ≥ 3 en posttest)
   │
   ├─ ev_obs ← ev_quad_from_data(n1, n2, x1, x2, a0, a1, a2)
   │  └─ Calcula e-valor observado del GRUPO REAL
   │
   ├─ ev_H ← simulate_evs_H_post(n1, n2, x1, x2, M=2000)
   │  └─ Simula 2000 datasets muestreando θ desde la POSTERIOR
   │      (dado que observamos x1, x2)
   │
   ├─ ev_A ← simulate_evs_A_post(n1, n2, x1, x2, M=2000)
   │  └─ Simula 2000 datasets desde la POSTERIOR completa
   │
   ├─ k_star ← find_kstar(ev_H, ev_A)
   │
   ├─ decision ← if(ev_obs ≤ k_star) "Reject H" else "Do not reject"
   │
   └─ Imprime fila de tabla con: grupo, n1, n2, x1, x2, ev_obs, k_star, decision

6. Salida: output/kstar_posterior_table.tex
```

**Por qué es diferente de Tabla A:** 
- Tabla A: ¿Cuál debería ser el cutoff ANTES de ver datos?
- Tabla B: Dado QUE vimos estos datos específicos, ¿cuál es el cutoff y qué decidimos?

---

## 4. PIPELINE C: GENERAR FIGURAS

**Archivo:** `make_figures.R`

**Objetivo:** Visualizar las curvas de error y superficies posteriores para el paper.

**Pasos:**

```
1. sourceCpp("BivBetaBinom.cpp")

2. FIGURAS PRIOR-BASED (réplica de figs 9-12 del paper):
   
   Para n ∈ {25, 50, 75, 100}:
   ├─ ev_H ← simulate_evs_H(n, n, M=1000)
   ├─ ev_A ← simulate_evs_A(n, n, M=1000)
   └─ plot: gráfico con 3 líneas
      ├─ α(k) = proporción de ev_H ≤ k
      ├─ β(k) = proporción de ev_A > k
      └─ α(k) + β(k) con línea vertical en k*
      
   → Figures/error_curves_prior_n25.png
   → Figures/error_curves_prior_n50.png
   → ... etc

3. FIGURAS POSTERIOR-BASED (4 grupos THKS):
   
   Para cada grupo (yy, yn, ny, nn):
   ├─ Extrae n1, n2, x1, x2 del TVSFP
   ├─ consts ← bb_constants(n1, n2, x1, x2, a0, a1, a2)
   ├─ ev_H ← simulate_evs_H_post(n1, n2, x1, x2, M=1000)
   ├─ ev_A ← simulate_evs_A_post(n1, n2, x1, x2, M=1000)
   └─ plot: curva de error posterior
      
      → Figures/error_curves_post_yy.png
      → Figures/error_curves_post_yn.png
      → ... etc

4. SUPERFICIES POSTERIORES (4 grupos):
   
   Para cada grupo:
   ├─ consts ← bb_constants(...)
   ├─ Evalúa densidad en grilla 80×80
   └─ plot 3D: superficie de θ1 vs θ2 vs densidad
      
      → Figures/posterior_surface_yy.png
      → Figures/posterior_surface_yn.png
      → ... etc

5. Tiempo total: ~10 minutos
```

---

## 5. PIPELINE D: VALIDACIÓN (Opcional, para desarrollo)

**Archivo:** `validate_ev_quad.R`

**Objetivo:** Verificar que la cuadratura Simpson 2D da MISMO resultado que MCMC (Stan).

**Pasos:**

```
1. sourceCpp("BivBetaBinom.cpp")
2. stan_model_obj ← stan_model("BBpost3.stan")

3. Para cada grupo THKS:
   ├─ Calcula ev mediante cuadratura Simpson 2D:
   │  ev_quad = ev_quad(consts, sup_H)
   │  tiempo: ~50ms
   │
   ├─ Calcula ev mediante MCMC/Stan:
   │  fit ← sampling(stan_model, data=..., iter=10000, chains=4)
   │  ev_mcmc = ev_FBST(posterior_samples)
   │  tiempo: ~60 segundos
   │
   └─ Compara: |ev_quad - ev_mcmc| < 0.001  ✓
      
      Si están cerca → cuadratura está validada
      Si no → hay un bug

4. Resultado: ✓ Validado (ambos dan ~0.0004 para grupo yy)
```

**Por qué NO se usa en el pipeline principal:**
- Es lento (Stan toma ~1 min por grupo)
- La cuadratura ya está validada
- Solo necesitas Stan si quieres debugging

---

## 6. PIPELINE E: SUPLEMENTO (Opcional)

**Archivo:** `paper_supplement.R`

**Objetivo:** Análisis adicional con prior informativa (α=10).

**Pasos:**

```
1. Corre TODO igual que Tabla B, PERO con:
   a0 = a1 = a2 = 10  (prior informativa)
   
2. Genera Tabla 3: resultados con prior informativa

3. Genera Tabla 4: prior informativa EN CONFLICTO con datos
   a0 = 90, a1 = a2 = 10  (prior muy concentrada en θ=0.1)

4. Genera figuras comparativas:
   ├─ Comparación FBST vs McNemar
   ├─ k* adaptativo vs tamaño muestral
   ├─ Sensibilidad a hiperparámetros
   └─ Ejemplo de no-rechazo

5. Tiempo: ~30 minutos
```

---

## 7. PIPELINE F: SIMULACIÓN ANTIGUA (NO SE USA AHORA)

**Archivo:** `THKS_run.R`

**Objetivo:** Análisis exploratorio tipo "demostración" (no para el paper).

**Pasos:**

```
1. sourceCpp("BivBetaBinom.cpp")

2. Prepara 4 grupos igual que antes

3. Para cada grupo:
   ├─ Calcula SUPERFICIE POSTERIOR:
   │  ├─ Grilla de θ1, θ2
   │  └─ Evalúa densBB_cpp en cada punto
   │  └─ plot 3D con persp()
   │
   ├─ Calcula MODO (por GA, algoritmo genético):
   │  GA_mod ← ga(type="real-valued", 
   │              fitness = function(v) densBB_cpp(v[1], v[2], consts))
   │  modo = GA_mod@solution
   │
   ├─ Calcula MEDIA (por integral2, cuadratura numérica de R):
   │  E[θ1] = integral2(function(t1, t2) t1 * densBB(t1, t2), 0, 1, 0, 1)
   │
   ├─ Calcula SUPREMO BAJO H (por GA 1D):
   │  sup_H = max densBB_H(θ) para θ ∈ [0, 1]
   │
   └─ [COMENTADO] Calcula e-valor por MCMC:
      ev ← stan(..., iter=10000)  ← Esto SÍ usa Stan, pero está comentado

4. Imprime resumen de resultados

5. Tiempo: ~2-5 minutos si no descomentas MCMC
```

**Por qué está aquí pero no se usa en el paper:**
- Es exploración, no análisis final
- Las figuras de "modo" y "media" no aparecen en el paper
- El MCMC comentado es solo para validación manual

---

## RESUMEN: ¿QUÉ CORRE PARA EL PAPER?

**INDISPENSABLE (estos SÍ se corren):**

1. `build_article_tables.R` 
   - Genera Tabla A y B
   - Toma: ~3 horas
   - Salida: 2 archivos .tex

2. `make_figures.R`
   - Genera 12 figuras
   - Toma: ~10 minutos
   - Salida: 12 archivos .png

**OPCIONAL (para validación/suplemento):**

3. `validate_ev_quad.R` 
   - Verifica cuadratura vs Stan
   - Toma: ~5 minutos
   - Solo si desconfías de la cuadratura

4. `paper_supplement.R`
   - Análisis con prior informativa
   - Toma: ~30 minutos
   - Para apéndice/material suplementario

**IGNORAR (antiguo/desarrollo):**

5. `THKS_run.R`
   - Era exploración inicial
   - Tiene MCMC comentado
   - No genera salida para el paper

---

## ¿DÓNDE ESTÁ EL CAMBIO DE MCMC A CUADRATURA?

**ANTES (metodología antigua):**
- Para calcular e-valor: `fit ← stan(..., iter=10000)` → caro
- Para cada simulación en Tabla A: 81 × 2000 = 162,000 llamadas a Stan

**AHORA (metodología actual):**
- Para calcular e-valor: `ev_quad(consts, sup_H)` → rápido (C++)
- Para cada simulación: evalúa `densBB_cpp` que es ~1000× más rápido

**El cambio está en:**
- `ev_quad()` línea 193-237 de `BivBetaBinom.cpp`
- `simulate_evs_H()` línea 329-344 usa `ev_quad`, no Stan
- `simulate_evs_A()` línea 348-362 usa `ev_quad`, no Stan

**Lo que no cambió:**
- La posterior sigue siendo la misma (ecuación \ref{posti})
- Las priors siguen siendo Olkin-Liu bivariadas
- FBST sigue siendo la misma teoría

Solo cambió **CÓMO se calcula**: de aproximado (MCMC) a exacto (cuadratura).

---

## ÚLTIMA COSA: FLUJO DEL USUARIO (TÚ)

Si quisieras reproducir el paper COMPLETO:

```
Paso 1: cd c:/Users/Tomas/Bayesian_pretest_posttest

Paso 2: Abre R y corre:
        source("build_article_tables.R")
        → Espera 3 horas
        → Verifica que existan:
          - output/kstar_prior_table.tex
          - output/kstar_posterior_table.tex

Paso 3: Corre:
        source("make_figures.R")
        → Espera 10 minutos
        → Verifica que existan figuras en Figures/

Paso 4: Edita article.tex y verifica:
        - Métodos completos
        - Conclusiones escritas
        - \input{output/kstar_prior_table.tex} funciona
        - \input{output/kstar_posterior_table.tex} funciona

Paso 5: Compila PDF:
        pdflatex article.tex
        bibtex article
        pdflatex article.tex
        pdflatex article.tex

Paso 6: ¡Listo!
```

¿Más claro así?
