/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Genetico.h
 * Author: PORTATIL
 *
 * Created on 30 de octubre de 2019, 16:02
 */

#ifndef GENETICO_H
#define GENETICO_H

#include <vector>
#include <string>
#include "Individuo.h"

using namespace std;

class Genetico {
public:
    Individuo genetico(string rutaLog, vector<vector<int>>& flu, 
            vector<vector<int>>& dis, bool sim);
private:
    // Variables
    vector<Individuo> poblacion;
    vector<Individuo> nuevaPoblacion;
    int evaluaciones = 0;
    vector<int> elite;
    vector<int> peor;
    Individuo mejor;
    
    // Funciones privadas
    void buscaMejor();
    void buscaElites();
    void sustituyeElites();
    void buscaPeores();
    Individuo torneo();
    vector<int> OX2(vector<int> padre, vector<int> madre, vector<int> marcados);
    vector<int> MOC(vector<int> padre, vector<int> madre, int marca);
    vector<int> mutacion(vector<int> sol, int pos);
};

#endif /* GENETICO_H */

