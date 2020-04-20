/* 
 * File:   Individuo.h
 * Author: sjm00010
 *
 * Created on 23 de octubre de 2019, 16:01
 */

#ifndef INDIVIDUO_H
#define INDIVIDUO_H

#include <vector>

using namespace std;

class Individuo {
public:
    Individuo() = default;
    Individuo(const Individuo& orig);
    virtual ~Individuo();
    bool IsEvaluado() const;
    vector<int> GetSolucion();
    int GetCoste() const;
    void SetSolucion(vector<int> solucion, int gene);
    void creaIndividuo(vector<int> solucion, int coste, int gene);
    void setCoste(int coste);
    int GetGeneracion() const;
private:
    vector<int> solucion;
    int coste;
    bool evaluado;
    int generacion; // Generacion en la que se creo
};

#endif /* INDIVIDUO_H */

