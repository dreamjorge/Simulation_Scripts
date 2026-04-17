# Documento técnico — análisis del ray tracing actual y plan de mejora

## 1. Qué problema resuelve este módulo

Este módulo no implementa ray tracing volumétrico general ni path tracing. Lo que hace es **trazado geométrico de rayos inducido por el gradiente de fase del campo óptico**. En términos físicos, está siguiendo la dirección local del vector de onda efectivo del haz a partir de la fase de `beam.opticalField(x, y, z)`.

La arquitectura actual toma un haz analítico o semianalítico, evalúa su campo complejo en posiciones cercanas y usa esa información para aproximar las pendientes del rayo. Después integra esas pendientes a lo largo de `z` con Euler o RK4. La idea base es correcta para óptica paraxial y visualización de frentes de onda, caústicas y divergencia, pero la implementación todavía mezcla una base física buena con varias decisiones numéricas que le bajan precisión y confiabilidad.

## 2. Cómo está calculado hoy

El flujo principal está en `src/propagation/rays/RayTracePropagator.m` y `src/propagation/rays/RayTracer.m`. `RayTracePropagator` arma un `RayBundle` sembrado sobre una grilla 2D y delega a `RayTracer.propagate(...)`. Si no se pasa `dz`, usa una discretización fija de `100` pasos en `z`. Eso simplifica mucho la API, pero no adapta el paso al comportamiento local del haz.

En `RayTracer.propagate`, cada iteración toma la última capa del bundle, calcula pendientes locales `sx` y `sy`, y actualiza `x` e `y`. Si el método es Euler, usa una sola evaluación del gradiente. Si el método es RK4, evalúa cuatro veces y combina los incrementos. La decisión de ofrecer RK4 es correcta porque, para un mismo `dz`, reduce bastante el error de trayectoria frente a Euler.

La parte delicada está en `RayTracer.calculateSlopes`. Ahí se calcula la fase con `unwrap(angle(field))`, luego se vuelve a evaluar el campo en `x + delta` y `y + delta`, y se estiman derivadas por diferencia hacia adelante. Finalmente se divide por `k = 2*pi/lambda`, lo cual implementa la relación paraxial entre gradiente de fase y pendiente del rayo. Conceptualmente está bien. Numéricamente, ahí vive el mayor cuello de botella de calidad.

## 3. Qué está bien diseñado

El acoplamiento con los haces está BIEN resuelto. `RayTracer` no conoce detalles internos de Gaussian, Laguerre o Hankel; solo exige que el beam implemente `opticalField(X, Y, z)`. Eso deja una puerta excelente para mejorar la calidad sin romper toda la arquitectura.

También está bien que exista una variante `HankelRayTracer`, porque los haces Hankel tienen semántica física especial. La idea de llevar `ht` por rayo para distinguir ramas de Hankel tiene sentido. El problema no es la intención, sino el criterio concreto usado para detectar cruce de eje.

## 4. Dónde están los problemas numéricos reales

El primer problema serio es el **gradiente de fase por diferencias forward con delta fijo**. `delta = 1e-7` está hardcodeado y no depende ni del waist, ni del grid, ni de la longitud de onda, ni de la escala espacial local. Eso significa que a veces el paso es demasiado grande para captar curvatura fina y otras veces es demasiado chico y entra en cancelación numérica. Además, la diferencia forward tiene error de truncamiento peor que una diferencia central.

El segundo problema es usar `unwrap(angle(field))` como paso previo a la derivada espacial. Eso funciona más o menos en casos suaves, pero se vuelve frágil cerca de singularidades de fase, vórtices, regiones de baja amplitud o cambios bruscos de rama. Cuando uno deriva una fase unwrapped sobre una malla 2D, el algoritmo de unwrap puede introducir decisiones locales que no siempre son coherentes con la derivada física que querés medir. Acá no alcanza con “que no dé NaN”. Necesitás consistencia física.

El tercer problema es estructural: `RayBundle` guarda `r` y `theta` al construir el objeto, pero `addStep()` no los vuelve a calcular. Eso deja el objeto en un estado parcialmente inconsistente. Si otro módulo asume que `bundle.r` o `bundle.theta` reflejan la última posición propagada, va a leer basura vieja. Este tipo de bug es insidioso porque no siempre rompe de inmediato; a veces solo degrada resultados.

El cuarto problema está en `HankelRayTracer`. Para decidir si un rayo cruzó el eje óptico, se usa el signo de `(x0*y1 - x1*y0) < 0`. Eso no demuestra que el segmento entre `(x0,y0)` y `(x1,y1)` pase por el origen. Ese valor solo mide orientación relativa entre vectores. O sea: la lógica puede flippear de `H^(2)` a `H^(1)` cuando no hubo cruce real. Físicamente, eso es demasiado peligroso para dejarlo así.

## 5. Qué dicen hoy los tests y qué NO dicen

Los tests actuales sirven como smoke tests. Verifican inicialización, que la propagación agregue pasos, que los valores sean finitos y que en algunos casos haya divergencia. Eso está bien como base mínima. El problema es que varios tests de edge cases setean `bundle.sx` y `bundle.sy`, pero `RayTracer.propagate()` no usa esas slopes iniciales; recalcula todo desde el campo en cada iteración. Entonces esos tests no están validando lo que dicen validar. Parecen cubrir backward slopes o zero slopes, pero en realidad no ejercitan esa lógica.

También falta validación analítica. Para GaussianBeam sí tenés una oportunidad de oro: su fase está formulada de manera explícita, así que podés comparar el gradiente numérico contra el gradiente analítico esperado. Sin esa capa, el módulo puede “correr” y seguir estando mal.

## 6. Técnica recomendada para mejorar el cálculo

La mejora más rentable no es cambiar toda la arquitectura. Es cambiar **cómo se calcula el gradiente de fase**.

La primera versión mejorada debería reemplazar la diferencia forward por **diferencia central**, con un `delta` escalado respecto de la geometría del problema. Eso solo ya reduce error y sesgo.

La segunda versión, que es la que realmente vale la pena si querés hacer esto bien, debería evitar derivar la fase unwrapped y trabajar directamente con el campo complejo. La formulación recomendada es extraer el gradiente de fase mediante una identidad del tipo:

\[
\nabla \phi = \mathrm{Im}\left(\frac{\nabla u}{u}\right)
\]

o, mejor expresada numéricamente en forma regularizada,

\[
\nabla \phi = \frac{\mathrm{Im}\left(\overline{u}\,\nabla u\right)}{|u|^2 + \epsilon}
\]

Eso te evita buena parte del dolor de `unwrap`, mejora la estabilidad cerca de discontinuidades de fase y hace que el cálculo sea más fiel al campo complejo que ya tiene el beam. Es una mejora técnica de verdad, no maquillaje.

## 7. Técnica recomendada para la evolución de la arquitectura

El paso siguiente es separar responsabilidades. Hoy `RayTracer` sabe demasiado sobre cómo obtener slopes. Lo mejor sería introducir una API opcional en los beams, algo como `phaseGradient(X, Y, z)` o `raySlopes(X, Y, z)`. Si un beam tiene fórmula analítica, la usa. Si no la tiene, se cae a un método numérico genérico.

Esa estrategia te da tres niveles. Gaussian puede tener derivada analítica. Hankel puede tener una implementación especializada. El resto puede usar el fallback numérico complejo. Eso te deja una arquitectura extensible y, lo más importante, verificable. CONCEPTOS > CODE, loco: cuando el diseño refleja la física, mantenerlo deja de ser una locura cósmica.

## 8. Orden correcto de implementación

Primero corregiría el núcleo numérico del gradiente. Después sanearía `RayBundle`. Después arreglaría el criterio de axis crossing en Hankel. Recién ahí metería tests físicos de exactitud. Si hacés los tests antes de arreglar esas tres cosas, vas a llenarte de ruido. Si arreglás primero y medís después, vas a poder demostrar mejora real.

## 9. Resultado esperado

Si hacés esos cambios, el ray tracing va a pasar de ser una herramienta útil de visualización con precisión incierta a una base numérica bastante más seria. No hace falta reescribir el sistema. Hace falta endurecer el núcleo matemático, limpiar el estado del bundle y dejar de validar solo que “no explota”.

## 10. Archivos y cambios concretos

### `src/propagation/rays/RayTracer.m`

Reemplazar `calculateSlopes()` actual por una estrategia central-difference o, mejor, por gradiente complejo. Hoy usa `delta` fijo, `phase = unwrap(angle(field))` y forward difference. Lo propuesto es agregar un helper privado `calculatePhaseGradientComplex(...)`, otro `resolveDelta(...)`, y dejar central difference como fallback si hiciera falta. También conviene agregar protección cuando `abs(field)` sea muy chico. El beneficio es más precisión, menos sensibilidad al unwrap y mejor estabilidad cerca de singularidades.

### `src/propagation/rays/RayBundle.m`

Corregir la inconsistencia de estado. La opción recomendada es eliminar `r` y `theta` como propiedades almacenadas y convertirlas en `Dependent`, calculadas desde `x` e `y`. La alternativa es recalcularlas dentro de `addStep()`. La primera opción tiene menos riesgo de bugs de sincronización; la segunda puede ser apenas más rápida en algunos usos, pero es bastante más frágil.

### `src/propagation/rays/HankelRayTracer.m`

Reemplazar la lógica de axis crossing. Hoy `crossed = (x0.*y1 - x1.*y0) < 0` no prueba un cruce real por el origen. Lo correcto es calcular la distancia mínima del segmento al origen y verificar además que la proyección caiga dentro del tramo. Con eso, el flip de `H^(2)` a `H^(1)` pasa a depender de geometría real y no de una heurística de orientación.

### `src/propagation/rays/RayTracePropagator.m`

Mejorar la política de `dz`. Hoy se usa `dz = (z_final - z0)/100` cuando no se especifica nada. Lo propuesto es introducir un modo fijo, uno automático y uno adaptativo mediante alguna propiedad tipo `StepPolicy = 'fixed'|'auto'|'adaptive'`. Con eso podés controlar error numérico sin romper compatibilidad.

### `src/beams/GaussianBeam.m`

Agregar una API analítica opcional, por ejemplo `phaseGradient(X, Y, z)` o `raySlopes(X, Y, z)`. Eso permite usar el beam gaussiano como referencia exacta tanto para producción como para testing. Sirve para medir calidad del gradiente y verificar que el ray tracing converge a una solución esperada.

### `src/beams/HankelLaguerre.m` y eventualmente `src/beams/HankelHermite.m`

Evaluar una API especializada de slopes o phase gradient si la formulación física se puede expresar sin volver inmantenible el código. Si eso se complica demasiado, conviene dejar la implementación específica para Gaussian y usar un fallback numérico mejorado para el resto. El tradeoff es claro: analítico da más exactitud, numérico genérico da menor complejidad.

### `tests/modern/test_RayTracing.m`

Agregar tests de exactitud física, no solo de sanidad. Casos concretos: comparar slope numérica vs slope analítica para Gaussian, comparar Euler vs RK4 con tolerancias, verificar simetría radial y estudiar convergencia al refinar `dz`. Un test que solo confirma que no hay `NaN` no demuestra que el cálculo esté bien.

### `tests/edge_cases/test_RayTracing_extreme.m`

Corregir tests que hoy seteaban `bundle.sx/sy` esperando que influyeran en la propagación. Como `RayTracer.propagate()` recalcula las pendientes desde el field, esos tests no representan el comportamiento real. Hay dos caminos: cambiar el código para aceptar slopes iniciales explícitas, o reescribir los tests para validar el gradiente del field real. La recomendación es la segunda, salvo que se quiera soportar manual slope seeding como feature.

### `tests/modern/test_HankelRayTracing.m`

Agregar tests que prueben de verdad el flip de `ht` cuando hay cruce real del eje, que no haya flip cuando no lo hay, que las trayectorias sean consistentes para distintos tipos de Hankel y que el criterio geométrico no sea sensible a casos marginales fuera de tolerancia.

## 11. Orden de trabajo recomendado

La secuencia más razonable es empezar por `RayTracer.m`, seguir con `RayBundle.m`, después reforzar `tests/modern/test_RayTracing.m`, luego corregir `HankelRayTracer.m` y recién entonces endurecer `tests/modern/test_HankelRayTracing.m`. Ese orden da mejora técnica real rápido y además deja evidencia medible de la mejora.

## 12. Recomendación concreta

Si querés hacerlo BIEN sin abrir demasiados frentes al mismo tiempo, el orden recomendado es: primero `RayTracer.m`, segundo `RayBundle.m`, tercero `tests/modern/test_RayTracing.m`, cuarto `HankelRayTracer.m` y quinto `tests/modern/test_HankelRayTracing.m`.

Ese orden te da el mayor retorno técnico con el menor caos. Primero resolvés el núcleo matemático, después saneás el modelo de datos y recién ahí validás con tests físicos más exigentes.
