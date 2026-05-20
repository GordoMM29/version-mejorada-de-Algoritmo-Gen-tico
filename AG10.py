import random
import blosum
import copy
import time

blosum62 = blosum.BLOSUM(62)
NFE = 0
start_time = time.time()

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    # Devolver como listas de caracteres individuales
    return [list(seq1), list(seq2), list(seq3)]

def crear_individuo():
    """Crea un individuo copiando las secuencias originales"""
    secuencias_originales = get_sequences()
    return [seq[:] for seq in secuencias_originales]

def crear_poblacion_inicial(n=10):
    """Crea población inicial con n individuos"""
    poblacion = []
    for _ in range(n):
        individuo = crear_individuo()
        poblacion.append(individuo)
    return poblacion

def igualar_longitud_secuencias(individuo, gap='-'):
    """Iguala la longitud de todas las secuencias del individuo"""
    if not individuo:
        return individuo
    
    max_len = max(len(seq) for seq in individuo)
    individuo_igualado = []
    for seq in individuo:
        if len(seq) < max_len:
            seq_igualada = seq + [gap] * (max_len - len(seq))
        else:
            seq_igualada = seq[:]
        individuo_igualado.append(seq_igualada)
    return individuo_igualado

def evaluar_individuo_blosum62(individuo):
    """Evalúa un individuo usando la matriz BLOSUM62"""
    global NFE
    NFE += 1
    score = 0
    n_seqs = len(individuo)
    
    if n_seqs < 2:
        return 0
    
    # Asegurar que todas las secuencias tengan la misma longitud
    seq_len = len(individuo[0])
    for seq in individuo:
        if len(seq) != seq_len:
            # Si no tienen la misma longitud, igualarlas
            individuo = igualar_longitud_secuencias(individuo)
            seq_len = len(individuo[0])
            break
    
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i+1, n_seqs):
                # Aquí es crítico: individuo[i] es una lista de caracteres
                # individuo[i][col] debe ser un caracter
                a = individuo[i][col]
                b = individuo[j][col]
                
                # Verificar que sean strings (caracteres)
                if isinstance(a, list):
                    print(f"ERROR ESTRUCTURAL: a es lista, debería ser caracter")
                    print(f"individuo[{i}][{col}] = {a}")
                    return -999999  # Penalización grande
                
                if a == '-' or b == '-':
                    score -= 4
                else:
                    try:
                        score += blosum62[a][b]
                    except KeyError:
                        # Caracter no encontrado (puede ser por minúsculas)
                        score -= 2
    return score

def seleccion_torneo(poblacion, scores, k=3):
    """Selección por torneo - asegura que k no sea mayor que la población"""
    poblacion_size = len(poblacion)
    if poblacion_size == 0:
        return []
    
    k = min(k, poblacion_size)
    seleccionados = []
    
    for _ in range(poblacion_size):
        # Seleccionar k individuos diferentes
        if poblacion_size >= k:
            indices = random.sample(range(poblacion_size), k)
        else:
            indices = random.choices(range(poblacion_size), k=k)
        
        mejor_idx = max(indices, key=lambda i: scores[i])
        seleccionados.append(copy.deepcopy(poblacion[mejor_idx]))
    
    return seleccionados

def cruzar_individuos_multi_punto(ind1, ind2, num_cortes=3):
    """Cruza con múltiples puntos de corte"""
    hijo1 = []
    hijo2 = []
    
    for seq1, seq2 in zip(ind1, ind2):
        # Obtener índices donde hay aminoácidos (no gaps)
        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        
        if len(aa_indices) < num_cortes + 1:
            # Secuencia muy corta, copiar directamente
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue
        
        # Seleccionar puntos de corte
        num_cortes_efectivo = min(num_cortes, len(aa_indices) - 1)
        if num_cortes_efectivo > 0:
            cortes = sorted(random.sample(aa_indices, num_cortes_efectivo))
        else:
            cortes = []
        
        def cruza_multipunto(seqA, seqB, cortes):
            # Extraer aminoácidos
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            
            # Construir nueva secuencia de aminoácidos
            nueva_aa = []
            inicio = 0
            usar_A = True
            
            for corte in cortes + [len(aaA)]:
                if usar_A:
                    nueva_aa.extend(aaA[inicio:corte])
                else:
                    nueva_aa.extend(aaB[inicio:corte])
                inicio = corte
                usar_A = not usar_A
            
            # Reconstruir con gaps
            resultado = []
            idx_aa = 0
            for caracter in seqA:
                if caracter == '-':
                    resultado.append('-')
                else:
                    if idx_aa < len(nueva_aa):
                        resultado.append(nueva_aa[idx_aa])
                        idx_aa += 1
                    else:
                        resultado.append(caracter)
            
            return resultado
        
        nueva_seq1 = cruza_multipunto(seq1, seq2, cortes)
        nueva_seq2 = cruza_multipunto(seq2, seq1, cortes)
        
        hijo1.append(nueva_seq1)
        hijo2.append(nueva_seq2)
    
    return hijo1, hijo2

def mutar_individuo_avanzado(individuo, tasa_mutacion=0.3):
    """Mutación avanzada con múltiples operaciones"""
    nuevo_individuo = []
    
    for secuencia in individuo:
        sec = secuencia[:]
        
        if random.random() < tasa_mutacion:
            tipo_mutacion = random.choice(['insertar', 'eliminar', 'mover'])
            
            if tipo_mutacion == 'insertar':
                num_gaps = random.randint(1, 2)
                for _ in range(num_gaps):
                    pos = random.randint(0, len(sec))
                    sec.insert(pos, '-')
                    
            elif tipo_mutacion == 'eliminar' and '-' in sec:
                gap_positions = [i for i, a in enumerate(sec) if a == '-']
                if gap_positions:
                    pos = random.choice(gap_positions)
                    sec.pop(pos)
                    
            elif tipo_mutacion == 'mover' and '-' in sec:
                gap_positions = [i for i, a in enumerate(sec) if a == '-']
                if gap_positions:
                    pos_orig = random.choice(gap_positions)
                    gap = sec.pop(pos_orig)
                    pos_nueva = random.randint(0, len(sec))
                    sec.insert(pos_nueva, gap)
        
        nuevo_individuo.append(sec)
    
    return nuevo_individuo

def cruzar_poblacion_avanzada(poblacion, scores, num_cortes=3, elitismo=2):
    """Cruza avanzada con elitismo y selección por torneo"""
    if len(poblacion) == 0:
        return []
    
    # Elitismo
    poblacion_con_scores = list(zip(poblacion, scores))
    poblacion_con_scores.sort(key=lambda x: x[1], reverse=True)
    
    elite_size = min(elitismo, len(poblacion))
    elite = [copy.deepcopy(ind) for ind, _ in poblacion_con_scores[:elite_size]]
    
    # Seleccionar padres
    padres = seleccion_torneo(poblacion, scores, k=min(3, len(poblacion)))
    
    nueva_poblacion = elite[:]
    
    # Generar hijos
    random.shuffle(padres)
    for i in range(0, len(padres)-1, 2):
        padre1 = padres[i]
        padre2 = padres[i+1]
        
        hijo1, hijo2 = cruzar_individuos_multi_punto(padre1, padre2, num_cortes)
        
        if random.random() < 0.7:
            hijo1 = mutar_individuo_avanzado(hijo1)
        if random.random() < 0.7:
            hijo2 = mutar_individuo_avanzado(hijo2)
        
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)
    
    # Completar si es necesario
    while len(nueva_poblacion) < len(poblacion):
        nueva_poblacion.append(copy.deepcopy(random.choice(padres)))
    
    return nueva_poblacion[:len(poblacion)]

def refrescar_poblacion(poblacion, scores, tasa_refresco=0.2):
    """Refresca la población reemplazando los peores individuos"""
    if len(poblacion) == 0:
        return poblacion
    
    n_refrescar = max(1, int(len(poblacion) * tasa_refresco))
    
    # Identificar peores individuos
    indices_ordenados = sorted(range(len(scores)), key=lambda i: scores[i])
    peores_indices = indices_ordenados[:n_refrescar]
    
    # Crear nuevos individuos
    for idx in peores_indices:
        nuevo = crear_individuo()
        nuevo = mutar_individuo_avanzado(nuevo, tasa_mutacion=0.5)
        poblacion[idx] = nuevo
    
    return poblacion

def eliminar_peores_adaptativo(poblacion, scores, generacion, max_generaciones, porcentaje_base=0.5):
    """Eliminación adaptativa de los peores individuos"""
    if len(poblacion) <= 2:
        return poblacion, scores
    
    factor = 1 - (generacion / max_generaciones)
    porcentaje = porcentaje_base * (0.5 + factor * 0.5)
    
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = max(2, int(len(poblacion) * porcentaje))
    
    ind_seleccionados = [copy.deepcopy(poblacion[i]) for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    
    return ind_seleccionados, scores_seleccionados

def validar_poblacion_sin_gaps(poblacion, originales):
    """Valida la integridad de las secuencias"""
    for i, individuo in enumerate(poblacion):
        for j, (seq, seq_orig) in enumerate(zip(individuo, originales)):
            seq_sin_gaps = [a for a in seq if a != '-']
            if seq_sin_gaps != seq_orig:
                print(f"Error en individuo {i}, secuencia {j}")
                return False
    return True

def obtener_best(scores, poblacion):
    """Obtiene el mejor individuo de la población"""
    idx_mejor = scores.index(max(scores))
    fitness_best = scores[idx_mejor]
    best = copy.deepcopy(poblacion[idx_mejor])
    return best, fitness_best

if __name__ == "__main__":
    print("Iniciando Algoritmo Genético Mejorado...")
    print("="*50)
    
    # Parámetros
    TAMANO_POBLACION = 20
    MAX_GENERACIONES = 100
    
    veryBest = None
    fitnessVeryBest = float('-inf')
    
    # Inicializar población
    poblacion = crear_poblacion_inicial(TAMANO_POBLACION)
    poblacion = [mutar_individuo_avanzado(ind, tasa_mutacion=0.2) for ind in poblacion]
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    
    # Evaluación inicial
    scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
    
    # Verificar estructura
    print("Verificando estructura inicial...")
    for i, ind in enumerate(poblacion[:3]):
        print(f"Individuo {i}: {len(ind)} secuencias, longitud primera secuencia: {len(ind[0])}")
        print(f"  Primer caracter: {ind[0][0] if ind[0] else 'None'}")
    
    sin_mejora = 0
    mejor_fitness_anterior = float('-inf')
    
    for generacion in range(MAX_GENERACIONES):
        # Cruza
        num_cortes = random.randint(2, 4)
        poblacion = cruzar_poblacion_avanzada(poblacion, scores, num_cortes=num_cortes, elitismo=3)
        
        # Mutación
        poblacion = [mutar_individuo_avanzado(ind, tasa_mutacion=0.2) for ind in poblacion]
        
        # Refresco periódico
        if generacion > 0 and generacion % 15 == 0:
            poblacion = refrescar_poblacion(poblacion, scores, tasa_refresco=0.15)
        
        # Igualar longitudes y evaluar
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum62(ind) for ind in poblacion]
        
        # Eliminación adaptativa
        poblacion, scores = eliminar_peores_adaptativo(poblacion, scores, generacion, MAX_GENERACIONES, porcentaje_base=0.6)
        
        # Obtener mejor
        best, fitness_best = obtener_best(scores, poblacion)
        
        # Verificar estancamiento
        if fitness_best > mejor_fitness_anterior + 0.1:
            mejor_fitness_anterior = fitness_best
            sin_mejora = 0
        else:
            sin_mejora += 1
        
        # Actualizar mejor global
        if fitness_best > fitnessVeryBest:
            veryBest = copy.deepcopy(best)
            fitnessVeryBest = fitness_best
            sin_mejora = 0
        
        # Mostrar progreso
        end_time = time.time()
        transcurrido = end_time - start_time
        print(f"Gen {generacion:3d}: fitness={fitnessVeryBest:8.2f}, NFE={NFE:5d}, "
              f"tiempo={transcurrido:5.2f}s, sin_mejora={sin_mejora:2d}")
        
        # Aumentar mutación si hay estancamiento
        if sin_mejora > 20:
            print("  -> Aumentando tasa de mutación por estancamiento")
            poblacion = [mutar_individuo_avanzado(ind, tasa_mutacion=0.5) for ind in poblacion]
            sin_mejora = 0
    
    # Resultados finales
    print("\n" + "="*50)
    print("=== MEJOR INDIVIDUO ENCONTRADO ===")
    for i, seq in enumerate(veryBest):
        seq_str = ''.join(seq)
        print(f"Secuencia {i+1}: {seq_str[:100]}...")
        print(f"  Longitud: {len(seq)}")
    
    print(f"\nFitness final: {fitnessVeryBest:.2f}")
    print(f"NFE total: {NFE}")
    print(f"Tiempo total: {time.time() - start_time:.2f} segundos")
    
    originales = get_sequences()
    print(f"Validación de integridad: {validar_poblacion_sin_gaps(poblacion, originales)}")