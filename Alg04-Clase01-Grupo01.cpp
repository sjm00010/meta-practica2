/* 
 * File:   Genetico.cpp
 * Author: PORTATIL
 * 
 * Created on 30 de octubre de 2019, 16:02
 */

#include "Alg04-Clase01-Grupo01.h"
#include "Parametros.h"
#include "FuncionesComunes.h"
#include "random.h"

#include <limits>
#include <iostream>

using namespace std;

Individuo Genetico::genetico(string rutaLog, vector<vector<int>>& flu, 
        vector<vector<int>>& dis, bool sim){
    
    // Inicio la generacion
    int generacion = 0;
    mejor.setCoste(numeric_limits<int>::max());
    elite.resize(stoi(parametros[TAM_ELITE]));
    peor.resize(stoi(parametros[TAM_ELITE]));

    // Redimensiono el vector
    poblacion.resize(stoi(parametros[TAM_POBLACION]));

    for(int i = 0; i < poblacion.size(); i++){
        vector<int> sol = creaSolucion(flu.size());
        int coste = calculaCoste(sol ,flu,dis,sim);
        poblacion[i].creaIndividuo(sol , coste, generacion);
        evaluaciones++;
    }

    // Registro en el log
    if(parametros[LOG] == "SI"){
        creaLog(rutaLog);
        registraLogGeneracion(rutaLog, generacion);
        registraLogIndividuos(rutaLog, poblacion);
    }
    
    buscaMejor();
    
    generacion++;
    
    nuevaPoblacion.resize(poblacion.size());

    while(evaluaciones < stoi(parametros[LIM_EVALUACIONES])){
        
        // Busco los élites de la generación de partida
        buscaElites();
        
        // Registro en el log la poblacion de partida
        if(parametros[LOG] == "SI"){
            if(generacion < 30 || generacion % 10 == 0){
                registraLogGeneracion(rutaLog, generacion);
                registraCadena(rutaLog, "Poblacion de partida");
                registraLogIndividuos(rutaLog,poblacion);
            }
        }

        // Creo la nueva generación, torneo
        for(int i = 0; i < poblacion.size(); i++){
            Individuo aux = torneo();
            nuevaPoblacion[i].creaIndividuo(aux.GetSolucion(), aux.GetCoste(), aux.GetGeneracion());
        }

        // Cruzo la nueva población
        if(parametros[TIPO_CRUCE] == "OX2"){
            for(int i = 0; i < nuevaPoblacion.size(); i += 2){
                if(Randfloat(0,1) <= stof(parametros[PROB_CRUCE])){
                    vector<int> marcados;
                    vector<int> padre = nuevaPoblacion[i].GetSolucion();
                    for(int j= 0; j < nuevaPoblacion.size(); j++){
                        if(Randfloat(0,flu.size()) < stoi(parametros[PROB_OX2])){
                            marcados.push_back(padre[j]);
                        }
                    }
                    
                    padre = OX2(nuevaPoblacion[i].GetSolucion(), nuevaPoblacion[i+1].GetSolucion(), marcados);
                    vector<int> madre = OX2(nuevaPoblacion[i+1].GetSolucion(), nuevaPoblacion[i].GetSolucion(), marcados);
                    
                    nuevaPoblacion[i].SetSolucion(padre, generacion);
                    nuevaPoblacion[i+1].SetSolucion(madre, generacion);
                }
            }
        }else if(parametros[TIPO_CRUCE] == "MOC"){
            for(int i = 0; i < nuevaPoblacion.size(); i += 2){
                if(Randfloat(0,1) <= stof(parametros[PROB_CRUCE])){
                    vector<int> padre = nuevaPoblacion[i].GetSolucion();
                    int marca = Randint(1,padre.size()-2);
                    
                    padre = MOC(nuevaPoblacion[i].GetSolucion(), nuevaPoblacion[i+1].GetSolucion(), marca);
                    vector<int> madre = MOC(nuevaPoblacion[i+1].GetSolucion(), nuevaPoblacion[i].GetSolucion(), marca);
                    
                    nuevaPoblacion[i].SetSolucion(padre, generacion);
                    nuevaPoblacion[i+1].SetSolucion(madre, generacion);
                }
            }
        }else{
            cerr << "Error en la seleccion del operador de cruce" << endl;
        }
        
        // Registro en el log la poblacion
        if(parametros[LOG] == "SI"){
            if(generacion < 30 || generacion % 10 == 0){ 
                registraCadena(rutaLog, "Poblacion despues del cruce");
                registraLogIntermedio(rutaLog,poblacion);
            }
        }

        
        // Muto la nueva población
        for(int i = 0; i < nuevaPoblacion.size(); i++){
            for(int j = 0; j < nuevaPoblacion[i].GetSolucion().size(); j++){
                if(Randfloat(0,1) < stof(parametros[PROB_MUTA])){
                    // Mutación
                    nuevaPoblacion[i].SetSolucion(mutacion(nuevaPoblacion[i].GetSolucion(), j), generacion);
                }
            }
        }
        
        if(parametros[LOG] == "SI"){
            if(generacion < 30 || generacion % 10 == 0){
                // Registro en el log la poblacion
                registraCadena(rutaLog, "Poblacion despues de la mutacion");
                registraLogIntermedio(rutaLog,poblacion);
            }
        }

        // Calculo los costes de los nuevos genes
        for(int i = 0; i < nuevaPoblacion.size(); i++){
            vector<int> prueba = nuevaPoblacion[i].GetSolucion();
            if(!nuevaPoblacion[i].IsEvaluado()){
                int coste = calculaCoste(nuevaPoblacion[i].GetSolucion(), flu, dis, sim);
                nuevaPoblacion[i].setCoste(coste);
                evaluaciones++;
            }
        }

        // Elites
        sustituyeElites();

        // Finalmente cambio la poblacion anterior por la nueva
        poblacion = nuevaPoblacion;
        generacion++;
        
        // Registro en el log la poblacion
        if(parametros[LOG] == "SI"){
            if(generacion < 30 || generacion % 10 == 0){
                registraCadena(rutaLog, "Se sustituye la poblacion nueva por la de partida");
            }
        }
        buscaMejor();

    }

    // Registro en el log la poblacion final
    registraLogGeneracion(rutaLog, generacion);
    registraCadena(rutaLog, "Poblacion de partida");
    registraLogIndividuos(rutaLog,poblacion);
            
    return mejor;
}

// Funciones locales complementarias

/**
 * Funcion que busca el mejor de la poblacion, devolviendolo y borrandolo 
 * del vector
 * @return Una copia del mejor individuo
 */
void Genetico::buscaMejor(){
    int coste = mejor.GetCoste();
    int sol = -1;

    for(int i = 0; i < poblacion.size(); i++){
        if(poblacion[i].GetCoste() < coste){
            coste = poblacion[i].GetCoste();
            sol = i;
        }
    }
    if (sol != -1){
        mejor = poblacion[sol];
    }
}

/**
 * Funcion que busca el mejor de la poblacion, devolviendolo y borrandolo 
 * del vector
 * @return Una copia del mejor individuo
 */
void Genetico::buscaElites(){
    for(int i = 0; i < elite.size(); i++){
        elite[i] = 0;
    }
    for(int i = 0; i < poblacion.size(); i++){
        for(int j = 0; j < elite.size(); j++){
            if(poblacion[i].GetCoste() < poblacion[elite[j]].GetCoste()){
                for(int k = j+1; k < elite.size();k++)
                    elite[k] = elite[k-1];
                elite[j] = i;
                j=elite.size(); // Fin del bucle j
            }
        }
    }
}

/**
 * Funcion que busca los peores de la poblacion
 */
void Genetico::buscaPeores(){
    for(int i = 0; i < peor.size(); i++){
        peor[i] = 0;
    }
    for(int i = 0; i < poblacion.size(); i++){
        for(int j = 0; j < peor.size(); j++){
            if(nuevaPoblacion[i].GetCoste() > nuevaPoblacion[peor[j]].GetCoste()){
                for(int k = j+1; k < peor.size();k++)
                    peor[k] = peor[k-1];
                peor[j] = i;
                j=peor.size(); // Fin del bucle j
            }
        }
    }
}



/**
 * Funcion que busca si los elites estan en la generacion actual
 */
void Genetico::sustituyeElites() {
    vector<bool> marca(elite.size());
    for(int i = 0; i < marca.size(); i++){
        marca[i] = false;
    }

    buscaPeores();
    
    // Busca si el individuo ya se encuentra en la nuevaPoblacion
    for(int i = 0; i < nuevaPoblacion.size(); i++){
        for(int j = 0; j < elite.size(); j++){
            if(nuevaPoblacion[i].GetSolucion() == poblacion[elite[j]].GetSolucion()){
                marca[j] = true;
                j = elite.size(); // Fin del bucle j
            }
        }
    }

    // Sustituye aquellos individuos que no estan
    for(int i = 0; i < elite.size(); i++){
        if(!marca[i]){
            nuevaPoblacion[peor[i]] = poblacion[elite[i]]; 
        } 
    }
}

/**
 * Funcion que realiza la seleccion del mejor de 2 individuos
 * @return El mejor de los individuos
 */
Individuo Genetico::torneo(){
    int ind1 = Randint(0, poblacion.size()-1);
    int ind2;
    while (ind1 == ind2){ind2 = Randint(0, poblacion.size()-1);}

    if(poblacion[ind1].GetCoste() < poblacion[ind2].GetCoste()){
        return poblacion[ind1];
    }

    return poblacion[ind2];
}

/**
 * Funcion que implementa el operador de cruce OX2
 * @param padre Primer gen a cruzar
 * @param madre Segundo gen a cruzar
 * @param marcados Vector de los elementos marcados para cruzar
 * @return Devuelve el hijo
 */
vector<int> Genetico::OX2(vector<int> padre, vector<int> madre, vector<int> marcados){

    vector<int> posiciones;
    vector<int> datos;
    int total = 1;
    for(int i = 0; i < padre.size() && total < (marcados.size()*2) ; i++){
        for(int j = 0; j < marcados.size(); j++){
            if(padre[i] == marcados[j]){
                posiciones.push_back(i);
                total++;
            }

            if(madre[i] == marcados[j]){
                datos.push_back(madre[i]);
                total++;
            }
        }
    }

    for (int i = 0; i < posiciones.size(); i++){
        padre[posiciones[i]] = datos[i];
    }
    return padre;
}

/**
 * Funcion que implementa el operador de cruce MOC
 * @param padre Primer gen a cruzar
 * @param madre Segundo gen a cruzar
 * @param marca Marca a partir de la que se produce el cruce
 * @return Devuelve el hijo
 */
vector<int> Genetico::MOC(vector<int> padre, vector<int> madre, int marca){
    int pos = marca;
    for(int i = 0; i < padre.size(); i++){
        // Busco si el elemento debe conservarse
        bool marcado = false;
        for(int j = 0; j < marca; j++){
            if(padre[i] == madre[j])
                marcado = true;
        }

        // Lo sustituyo en caso de que no deba conservarse
        if(!marcado){
            padre[i] = madre[pos];
            pos++;
        }
    }
    return padre;
}


/**
 * Funcion que realiza la mutacion de los genes
 * @param sol Vector que hay que mutar
 * @param pos Posicion de la mutacion
 * @return Devuelve el vector mutado
 */
vector<int> Genetico::mutacion(vector<int> sol, int pos){
    int pos2, pos3;
    do{
        pos2 = Randint(0, sol.size()-1);
        pos3 = Randint(0, sol.size()-1);
    }while(pos != pos2 && pos2 != pos3 && pos != pos3 && pos2 < sol.size() 
            && pos3 < sol.size());

    swap(sol[pos],sol[pos3]);
    swap(sol[pos2],sol[pos3]);

    return sol;
}


