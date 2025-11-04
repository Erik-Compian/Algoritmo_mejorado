# Algoritmo_mejorado
Versión mejorada de un algoritmo genético con elitismo, optimización de mutación y comparación con versión original (PROYECTO FINAL)

Proyecto Final: Algoritmo Genético Mejorado

Este proyecto implementa una versión mejorada de un algoritmo genético para la alineación de secuencias biológicas utilizando la matriz BLOSUM62.
El objetivo fue optimizar el proceso evolutivo mediante la incorporación de nuevas estrategias de selección, cruza y mutación.

Características principales

Versión original: Algoritmo básico de alineación con cruza doble y mutación simple.

Versión mejorada:

Implementación de elitismo (los mejores individuos se preservan).

Mutación más eficiente con control de probabilidad.

Eliminación más estricta de los individuos menos aptos.

Comparación visual del fitness entre la versión original y la mejorada.

Resultados

La siguiente gráfica muestra la comparación del fitness obtenido por ambas versiones:

Como se puede observar, la versión mejorada con elitismo alcanza un mayor valor de fitness y una evolución más estable, lo que demuestra la eficacia de las mejoras aplicadas.

Validación de integridad

El código mantiene la validación de integridad para asegurar que las secuencias sin gaps no se alteren a lo largo del proceso evolutivo:

print("Validación de integridad:", validar_poblacion_sin_gaps(poblacion, originales))


Esto garantiza que la optimización no modifica los datos biológicos originales.

Conclusiones

El elitismo mejora significativamente la convergencia del algoritmo.

La mutación controlada permite mantener diversidad sin alterar la estructura base.

La selección optimizada acelera la evolución hacia soluciones de mejor calidad.

Este proyecto demuestra cómo pequeñas modificaciones en los operadores genéticos pueden generar mejoras sustanciales en el rendimiento global del algoritmo.

Autor

Nombre: Erik Alejandro Compian Ovalle
Materia: Análisis y Modelación de Sistemas / Proyecto final
Profesor: Dr. Ernesto Rios Willars.
