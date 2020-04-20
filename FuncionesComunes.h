/* 
 * File:   FuncionesComunes.h
 * Author: PORTATIL
 *
 * Created on 1 de noviembre de 2019, 21:39
 */

#ifndef FUNCIONESCOMUNES_H
#define FUNCIONESCOMUNES_H

#include "Individuo.h"

using namespace std;
    // Carga Datos
    void cargaParametros();
    void cargaDatos(string &archivo, vector<vector<int>> &flu,
            vector<vector<int>> &dis, bool &sim);
    
    // Registros .log
    void creaLog(string nombreAr);
    void registraCadena(string log, string elec);
    void registraTiempo(string log, double t, int semilla);
    void registraLogGeneracion(string log, int generacion);
    void registraLogIntermedio(string log, vector<Individuo> poblacion);
    void registraLogIndividuos(string log, vector<Individuo> poblacion);
    void registraLogSolucion(string log, Individuo i);
    
    // Otros
    vector<int> creaSolucion(int tam);
    int calculaCoste (vector<int> sol, vector<vector<int>>& flu, 
            vector<vector<int>>& dis, bool sim);
    int calculaSemilla(int prueba);
    void mostrarResultado( vector<int> v, double t, vector<vector<int>>& flu,
            vector<vector<int>>& dis, const int coste);

#endif /* FUNCIONESCOMUNES_H */

