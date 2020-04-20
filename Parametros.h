/*
 * File:   Parametros.h
 * Author: sjm00010
 *
 * Created on 27 de septiembre de 2019, 10:58
 */

#ifndef PARAMETROS_H
#define PARAMETROS_H

#include <string>

using namespace std;

// Parametros necesarios que se necesitan
extern string rutaParam;
extern int numParam;
extern vector<string> parametros;

// Parametros a cargar del archivo parametros.txt
enum valor {
    CARPETA_DATOS = 0,
    CARPETA_LOG = 1,
    NOMBRE_ARCHIVO = 2,
    DNI = 3,
    NUM_PRUEBAS = 4,
    TAM_POBLACION = 5,
    LIM_EVALUACIONES = 6,
    PROB_CRUCE = 7,
    PROB_MUTA = 8,
    PROB_OX2 = 9,
    TAM_ELITE = 10,
    TIPO_CRUCE = 11,
    LOG = 12,
};

#endif /* PARAMETROS_H */

