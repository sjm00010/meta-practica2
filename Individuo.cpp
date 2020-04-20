/* 
 * File:   Individuo.cpp
 * Author: PORTATIL
 * 
 * Created on 23 de octubre de 2019, 16:01
 */

#include "Individuo.h"

Individuo::Individuo(const Individuo& orig) {
    this->evaluado = orig.evaluado;
    this->solucion = orig.solucion;
    this->coste = orig.coste;
    this->generacion = orig.generacion;
}

Individuo::~Individuo() {
}

void Individuo::creaIndividuo(vector<int> solucion, int coste, int gene) {
    this->solucion = solucion;
    this->coste = coste;
    this->generacion = gene;
    this->evaluado = true;
}

void Individuo::setCoste(int coste) {
    this->coste = coste;
    this->evaluado = true;
}

int Individuo::GetGeneracion() const {
    return generacion;
}

bool Individuo::IsEvaluado() const {
    return evaluado;
}

vector<int> Individuo::GetSolucion() {
    return solucion;
}

int Individuo::GetCoste() const {
    return coste;
}

void Individuo::SetSolucion(vector<int> solucion, int gene) {
    this->solucion = solucion;
    this->evaluado = false;
    this->generacion = gene;
}



