import random
import pandas as pd
import blosum
import copy
import time
import matplotlib.pyplot as plt

# CONFIGURACIÓN GLOBAL

blosum62 = blosum.BLOSUM(62)
random.seed(42)  # reproducibilidad


# FUNCIONES BASE

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]


def crear_individuo():
    return get_sequences()


def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    poblacion = [[row[:] for row in individuo_base] for _ in range(n)]
    return poblacion


def mutar_poblacion_v2(poblacion, num_gaps=1):
    poblacion_mutada = []
    for individuo in poblacion:
        nuevo_individuo = []
        for fila in individuo:
            fila_mutada = fila[:]
            posiciones = set()
            for _ in range(num_gaps):
                pos = random.randint(0, len(fila_mutada))
                while pos in posiciones:
                    pos = random.randint(0, len(fila_mutada))
                posiciones.add(pos)
                fila_mutada.insert(pos, '-')
            nuevo_individuo.append(fila_mutada)
        poblacion_mutada.append(nuevo_individuo)
    return poblacion_mutada


def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    individuo_igualado = [fila + [gap] * (max_len - len(fila)) for fila in individuo]
    return individuo_igualado


def evaluar_individuo_blosum62(individuo):
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i + 1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += blosum62[a][b]
    return score


def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = int(len(poblacion) * porcentaje)
    ind_seleccionados = [poblacion[i] for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    return ind_seleccionados, scores_seleccionados


def mutar_individuo(individuo, n_gaps, p):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < p:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec))
                while pos in posiciones:
                    pos = random.randint(0, len(sec))
                posiciones.add(pos)
                sec.insert(pos, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo

# FUNCIÓN CRUZAR CORREGIDA

def cruzar_individuos_doble_punto(ind1, ind2):
    hijo1, hijo2 = [], []
    for seq1, seq2 in zip(ind1, ind2):
        # Asegurar que sean listas
        seq1 = list(seq1)
        seq2 = list(seq2)

        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        if len(aa_indices) < 6:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue

        p1, p2 = sorted(random.sample(aa_indices, 2))

        def cruza(seqA, seqB):
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            nueva = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
            resultado = []
            idx = 0
            for a in seqA:
                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva[idx])
                    idx += 1
            return resultado

        nueva_seq1 = cruza(seq1, seq2)
        nueva_seq2 = cruza(seq2, seq1)

        # Asegurar que sean listas antes de mutar
        nueva_seq1 = list(nueva_seq1)
        nueva_seq2 = list(nueva_seq2)

        # Corregido
        hijo1.append(mutar_individuo([nueva_seq1], 1, 0.8)[0])
        hijo2.append(mutar_individuo([nueva_seq2], 1, 0.8)[0])

    return hijo1, hijo2


def cruzar_poblacion_doble_punto(poblacion):
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    parejas = [(indices[i], indices[i + 1]) for i in range(0, n - 1, 2)]
    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))

    for idx1, idx2 in parejas:
        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]
        hijo1, hijo2 = cruzar_individuos_doble_punto(padre1, padre2)
        nueva_poblacion.extend([copy.deepcopy(padre1), copy.deepcopy(padre2), hijo1, hijo2])
    return nueva_poblacion[:2 * n]


def validar_poblacion_sin_gaps(poblacion, originales):
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True


def obtener_best(scores, poblacion):
    idx_mejor = scores.index(max(scores))
    return copy.deepcopy(poblacion[idx_mejor]), scores[idx_mejor]

# FUNCIÓN PRINCIPAL DE EJECUCIÓN

def ejecutar_algoritmo(usando_elitismo=False, generaciones=100):
    poblacion = crear_poblacion_inicial(10)
    poblacion = mutar_poblacion_v2(poblacion, num_gaps=1)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores)

    fitness_por_generacion = []
    originales = get_sequences()
    best_global = None
    best_fitness = float('-inf')

    for gen in range(generaciones):
        poblacion = cruzar_poblacion_doble_punto(poblacion)
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        poblacion, scores = eliminar_peores(poblacion, scores)

        best, fit = obtener_best(scores, poblacion)

        # Aplicar elitismo si está activado
        if usando_elitismo:
            poblacion[0] = best
            scores[0] = fit

        if fit > best_fitness:
            best_global = best
            best_fitness = fit

        fitness_por_generacion.append(best_fitness)

    print("Validación de integridad:", validar_poblacion_sin_gaps(poblacion, originales))
    return fitness_por_generacion, best_fitness

# EJECUCIÓN Y GRÁFICA COMPARATIV

if __name__ == "__main__":
    print("Ejecutando algoritmo original...")
    fitness_original, final_fit_orig = ejecutar_algoritmo(usando_elitismo=False)

    print("\nEjecutando algoritmo mejorado (con elitismo)...")
    fitness_mejorado, final_fit_mejorado = ejecutar_algoritmo(usando_elitismo=True)

    # Mostrar resultados
    print(f"\nFitness final - Original: {final_fit_orig}")
    print(f"Fitness final - Mejorado: {final_fit_mejorado}")

    # Graficar comparación
    plt.figure(figsize=(10, 6))
    plt.plot(fitness_original, label="Algoritmo original", color='blue')
    plt.plot(fitness_mejorado, label="Algoritmo mejorado (elitismo)", color='red')
    plt.title("Comparación del Fitness: Original vs Mejorado")
    plt.xlabel("Generaciones")
    plt.ylabel("Fitness (más alto es mejor)")
    plt.legend()
    plt.grid(True)

    # Guardar imagen para GitHub
    plt.savefig("comparacion_fitness.png", dpi=300)
    plt.show()

