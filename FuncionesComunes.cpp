/* 
 * File:   FuncionesComunes.cpp
 * Author: PORTATIL
 * 
 * Created on 1 de noviembre de 2019, 21:39
 */

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include "FuncionesComunes.h"
#include "Parametros.h"
#include "random.h"

using namespace std;

    /**
     * Función para cargar los parametros.
     * @param nombreArchivo Ruta del archivo con los parametros.
     */
    void cargaParametros() {
        ifstream fe;
        try {
            fe.open(rutaParam);
            string linea;
            while (!fe.eof()) {
                fe.ignore(256, ' ');
                fe.ignore(256, ' ');
                fe >> linea;
                parametros.push_back(linea.c_str());
            }
        } catch (ifstream::failure e){
            cerr << "Error en la lectura del archivo de parametros, verifique rutaParametros" << endl;
        }
    }
    
    /**
     * Función que comprueba tanto a simetria como la correcta atribución de las matrices
     * @param flu Matriz de flujos
     * @param dis Matriz de distancias
     * @return True si es simétrica, False si no lo es.
     */
    bool comprueba(vector<vector<int>> &flu,
            vector<vector<int>> &dis) {
        bool cambio = false;
        bool sime = true;
        for (int i = 0; i < dis.size(); i++) {
            for (int j = i + 1; j < dis.size(); j++) {
                if (dis[i][j] == 0 || dis[j][i] == 0) {
                    cambio = true;
                }
                if (dis[i][j] != dis[j][i]) {
                    sime = false;
                    // Finalizo el bucle
                    i = dis.size();
                    j = dis.size();
                }
            }
        }

        if (cambio) {
            //Cambio las matrices
            vector<vector<int>> temp(flu);
            flu = dis;
            dis = temp;
        }
        return sime;
    }
    
    /**
     * Función que carga los datos en las matrices.
     * @param archivo Ruta del archivo con los datos.
     * @param flu Matriz de flujos.
     * @param dis Matriz de distancias.
     * @param sim Variable que almacena si las matrices son simetricas.
     */
    void cargaDatos(string &archivo, vector<vector<int>> &flu,
            vector<vector<int>> &dis, bool &sim) {
        ifstream fe(archivo); //Creo un flujo de entrada

        if (fe.good()) {
            int tam;
            fe >> tam; // Recojo el tamaño de las matrices

            // Con el tamaño, redimensiono las matrices para adaptarlas al problema.
            flu.resize(tam);
            dis.resize(tam);
            for (int i = 0; i < tam; i++) {
                flu[i].resize(tam);
                dis[i].resize(tam);
            }

            // Leo los datos del fichero y los asigno a la matriz
            for (int i = 0; i < tam; i++) {
                for (int j = 0; j < tam; j++) {
                    fe >> dis[i][j];
                }
            }

            fe >> ws; //getline(cin, ws);

            for (int i = 0; i < tam; i++) {
                for (int j = 0; j < tam; j++) {
                    fe >> flu[i][j];
                }
            }
        }
        fe.close();

        // Compruebo las matrices
        sim = comprueba(flu, dis);
    }
    
    /**
     * Función que crea el archivo .log
     * @param nombre Nombre del archivo
     */
    void creaLog(string nombreAr) {
        // Crea un fichero de salida
        ofstream fs(nombreAr);
        fs.close();
    }
    
    /**
     * Función para añadir una cadena de texto al log
     * @param elec Cadena a añadir
     */
    void registraCadena(string log, string elec){
        ofstream fs;
        fs.open(log, std::ios_base::app);
        fs << "\n" << elec << "\n";
    }
    
    /**
     * Funcion que registra la generacion actual
     * @param log Ruta del archivo
     * @param generacion Generacion
     */
    void registraLogGeneracion(string log, int generacion){
        ofstream fs;
        fs.open(log, std::ios_base::app);
        fs << "----------------------------------------------------\n";
        fs << "              GENERACION " << generacion << "\n";
        fs << "----------------------------------------------------\n";
        fs.close();
    }
    
    /**
     * Funcion que registra los individuos tras el cruce y la mutacion
     * @param log Ruta del archivo
     * @param poblacion Poblacion a registrar
     */
    void registraLogIntermedio(string log, vector<Individuo> poblacion){
        ofstream fs;
        fs.open(log, std::ios_base::app);
        for(int i = 0; i < poblacion.size(); i++){
            fs << " Individuo " << i+1 << " :\n       ";
            vector<int> v = poblacion[i].GetSolucion();
            for (int j = 0; j < v.size(); j++) {
                fs << "(" << j + 1 << ")->" << v[j]+1 << "   ";
                if (j % 9 == 0 && j > 0) {
                    fs << "\n       ";
                }
            }
            fs << "\n";
        }
        fs.close();
    }
    
    /**
     * Funcion que registra todos los datos de la poblacion
     * @param log Ruta del archivo
     * @param poblacion Poblacion a registrar
     */
    void registraLogIndividuos(string log, vector<Individuo> poblacion){
        ofstream fs;
        fs.open(log, std::ios_base::app);
        for(int i = 0; i < poblacion.size(); i++){
            fs << " Individuo " << i+1 << " :\n       ";
            vector<int> v = poblacion[i].GetSolucion();
            for (int j = 0; j < v.size(); j++) {
                fs << "(" << j + 1 << ")->" << v[j]+1 << "   ";
                if (j % 9 == 0 && j > 0) {
                    fs << "\n       ";
                }
            }
            fs << "\n";
            fs << "        Coste : " << poblacion[i].GetCoste() << "\n";
            fs << "        Generacion de creacion : " << poblacion[i].GetGeneracion() << "\n";
        }
        fs.close();
    }
    
    void registraLogSolucion(string log, Individuo i){
        ofstream fs;
        fs.open(log, std::ios_base::app);
        fs << "----------------------------------------------------\n";
            fs << " Mejor Individuo :\n       ";
            vector<int> v = i.GetSolucion();
            for (int j = 0; j < v.size(); j++) {
                fs << "(" << j + 1 << ")->" << v[j]+1 << "   ";
                if (j % 9 == 0 && j > 0) {
                    fs << "\n       ";
                }
            }
            fs << "\n";
            fs << "        Coste : " << i.GetCoste() << "\n";
            fs << "        Generacion de creacion : " << i.GetGeneracion() << "\n";
        fs.close();
    }
    
    /**
     * Función para registrar el tiempo del algoritmo.
     * @param log Ruta del archivo .log
     * @param t Tiempo
     */
    void registraTiempo(string log, double t, int semilla) {
        ofstream fs;
        fs.open(log, std::ios_base::app);
        fs << "----------------------------------------------------\n";
        fs << "Tiempo :  " << t << " ms \n";
        fs << "Semilla : " << semilla << "\n";
    }
    
    /**
     * Función para crear una solución aleatoria de partida
     * @param tam Tamaño del vector solución
     * @return Devuelve el vector solucion creado.
     */
    vector<int> creaSolucion(int tam){
        vector<int> solucion(tam);
        for(int i = 0; i < tam; i++){
            solucion[i]=i;
        }  

        for(int i = 0; i < tam; i++){
            // La primera generacion de numeros aleatorios sale generan el mismo
            // número, por eso se deben de ejecutar un minimo de 2 veces el swap
            swap(solucion[Randint(0,tam-1)],solucion[Randint(0,tam-1)]);
        }
        
        return solucion;
    }
    
    /**
     * Función que calcula el coste de una solción.
     * @param sol Vector con la solución.
     * @param flu Matriz de flujos.
     * @param dis Matriz de distancias.
     * @return Devuelve el coste calculado.
     */
    int calculaCoste (vector<int> sol, vector<vector<int>>& flu,
            vector<vector<int>>& dis, bool sim){
        int coste = 0;
        if(sim){
            for (int i = 0; i < flu.size(); i++){
                for(int j = i+1; j < flu.size(); j++){
                    coste += flu[i][j]*dis[sol[i]][sol[j]];
                }
            }
            coste *= 2;
        }else{
            for (int i = 0; i < flu.size(); i++){
                for(int j = 0; j < flu.size(); j++){
                    if(i != j){
                        coste += flu[i][j]*dis[sol[i]][sol[j]];
                    }
                }
            }
        }
        
        return coste;
    }
    
    /**
     * Función que cambia la semilla automaticamente
     * @param prueba Numero de las 5 pruebas por el se se esta ejecutando
     * @return Nueva semilla
     */
    int calculaSemilla(int prueba){
        int dni = stoi(parametros[DNI], nullptr, 10);
        int d = dni > 0 ? (int) log10 ((double) dni) + 1 : 1;
		d = d - --prueba;
        d = pow(10,d);
        return (int)(dni%d*pow(10,prueba)+(int)dni/d);
    }
    
    /**
     * Función que muestra por consola la solción calculada.
     * @param v Vector con la solución.
     * @param t Tiempo que tarda el algoritmo.
     * @param flu Matriz de flujos
     * @param dis Matriz de distancias
     */
    void mostrarResultado( vector<int> v, double t, vector<vector<int>>& flu, 
            vector<vector<int>>& dis, const int coste){
        for(int i = 0; i < v.size(); i++){
            printf(" (%2d) -%3d  ",i+1, v[i]+1);
            if (i % 5 == 0 && i > 0 && i != v.size()-1){
                cout << "\n";
            }
        } 
        cout << "\n\n";
        cout << "   El coste de la solucion es : " << coste << "\n";
        cout << "   Tiempo empleado : " << t << " ms\n";
    }