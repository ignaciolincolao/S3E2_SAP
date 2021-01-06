#include <iostream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdint>
#include <stdio.h>
#include <sstream>
#include <random>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <string>

int n_students, n_colegios;

///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
struct Info_colegio{
    double latitude = 0.0;
    double longitude = 0.0;
    int num_alu = 0;
    int rbd = 0;
    int prioritario = 0;
};
///////////////////////////////////////////////////
/// Estructura de alumnos
///////////////////////////////////////////////////
struct Info_alu{
    int rbd = 0;
    int sep = 0;
    double latitude = 0.0;
    double longitude = 0.0;
};
///////////////////////////////////////////////////
/// Funciones generales
///////////////////////////////////////////////////
double calCosto(int currentSolution[], double **distMat, const double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]);
double meanDist(const int currentSolution[], double  **distMat);
double S(const int currentSolution[],const int alumnosSep[], int totalVuln);
double costCupo(const int currentSolution[],const int cupoArray[]);
void newSolution(int currentSolution[],const int previousSolution[]);
int acepta(double costPrevious, double costCurrent);
double p(double costPrevious,double costCurrent);
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]);
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat);


///////////////////////////////////////////////////
/// Parametros de configuración Default
///////////////////////////////////////////////////


double coolingRate = 0.98; // Tasa de enfriamiento
double temp = 10000000; // Temperatura Inicial
double min_temp =0.00000000009;// 0.00000009; // Minima temperatura que puede llegar
double alpha1 = 15; // Alpha de distancia
double alpha2 = 30; // Alpha de segregación
double alpha3 = 25; // Alpha de costocupo
double max_temp = pow(10,300);
double k_recalentamiento = 0.994;
int seed=rand();

std::string ruta_save = "./save/"; // Ruta para guardar los archivos


///////////////////////////////////////////////////
/// Variables globales.
///////////////////////////////////////////////////

double alpha[3]={alpha1,alpha2,alpha3}; // Valores del alpha con orden Distancia, Segregación, Costo Cupo
std::random_device rd;
std::mt19937 mt(rd());
std::uniform_int_distribution<int> dist(0,0);
std::uniform_int_distribution<int> dist2(0,0);
double max_dist=0.0;
double min_dist=0.0;
double init_dist=0.0;

///////////////////////////////////////////////////
///////////////////////////////////////////////////
int main(int argc, char *argv[]) {
    time_t hora_actual;
    struct tm * time_info;
    char timestr[20];
    int x=0,z=0;
    int totalVuln=0;

    ///////////////////////////////////////////////////
    /// Genera archivo de almacenamiento de datos
    ///////////////////////////////////////////////////
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);

    std::string prefijo_save = std::string(timestr);
    if (argc>1) {
        alpha1 = std::stod(argv[1]); // Alpha de distanciathreadIdx.x
        alpha2 = std::stod(argv[2]); // Alpha de segregación
        alpha3 = std::stod(argv[3]); // Alpha de costocupo
        alpha[0]=alpha1;
        alpha[1]=alpha2;
        alpha[2]=alpha3;
        coolingRate = std::stod(argv[4]); // Tasa de enfriamiento
        k_recalentamiento = std::stod(argv[5]);
        temp = std::stod(argv[6]); // Temperatura inicial
        min_temp = std::stod(argv[7]); // Minima temperatura que puede llegar
        ruta_save = argv[8];
        prefijo_save = argv[9];
        max_temp= pow(10,300);
        seed= std::stoi(argv[10]);
    }
    mt.seed(seed);

    std::ofstream info;
    std::string infotxt = "./save/" + std::string(timestr) +"-info.txt"; // concatenar
    info.open(infotxt);


    Info_colegio *ptr_colegios;
    Info_alu *ptr_students;

    ///////////////////////////////////////////////////
    /// Datos colegios
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura Info_colegio
    ///////////////////////////////////////////////////
    std::string line_colegios;
    std::ifstream info_school("colegios_utm.txt"); // concatenar
    std::getline(info_school, line_colegios);
    n_colegios = std::stoi(line_colegios);
    Info_colegio colegios[n_colegios];
    ptr_colegios = &colegios[0];
    while (std::getline(info_school, line_colegios)) {
        std::stringstream linestream(line_colegios);
        std::string data;
        std::getline(linestream, data, ',');
        ptr_colegios->rbd = std::stoi(data);
        std::getline(linestream, data, ',');
        ptr_colegios->latitude = std::stod(data);
        std::getline(linestream, data, ',');
        ptr_colegios->longitude = std::stod(data);
        std::getline(linestream, data, ',');
        ptr_colegios->num_alu = std::stoi(data);
        std::getline(linestream, data, ',');
        ptr_colegios->prioritario = std::stoi(data);
        ptr_colegios++;
    }

    ptr_colegios = &colegios[0]; // vuelve el puntero al inicio
    info_school.close();
    ///////////////////////////////////////////////////
    /// Datos Alumnos
    /// Lee el archivo linea por linea y luego lo agrega al arreglo de estructura info_student
    ///////////////////////////////////////////////////
    std::string line_student;
    std::ifstream info_student("alumnos_utm.txt"); // concatenar


    std::getline(info_student, line_student);
    n_students = std::stoi(line_student);
    Info_alu students[n_students];
    ptr_students = &students[0];
    while (std::getline(info_student, line_student)) {
        std::stringstream linestream(line_student);
        std::string data;
        std::getline(linestream, data, ',');
        ptr_students->rbd = std::stoi(data);
        std::getline(linestream, data, ',');
        ptr_students->latitude = std::stod(data);
        std::getline(linestream, data, ',');
        ptr_students->longitude = std::stod(data);
        std::getline(linestream, data, ',');
        ptr_students->sep = std::stoi(data);
        if (ptr_students->sep == 1) {
            totalVuln++;
        }
        ptr_students++;
    }
    ptr_students = &students[0]; // vuelve el puntero al inicio
    info_student.close();

    ///////////////////////////////////////////////////
    /// Inicializa Variables y arreglos
    ///////////////////////////////////////////////////


    int *previousSolution= nullptr;
    int *bestSolution= nullptr;
    int *currentSolution=nullptr;
    double **distMat=nullptr;
    int *cupoArray=nullptr;
    int *alumnosSep=nullptr;
    int count=0;

    previousSolution = (int *)malloc(sizeof(int)*n_students);
    bestSolution=(int *)malloc(sizeof(int)*n_students);
    currentSolution=(int *)malloc(sizeof(int)*n_students);
    distMat=(double **)malloc(sizeof(double)*n_students);
    cupoArray=(int *)malloc(sizeof(int)*n_colegios);
    alumnosSep = (int *)malloc( sizeof(int)*n_students);

    ///////////////////////////////////////////////////
    /// Se crear un arreglo donde el el valor es la posición del estudiante sep
    ///////////////////////////////////////////////////

    for(x=0; x < n_students; x++) {
        distMat[ x ]=(double *)malloc(sizeof(double)*n_colegios);
        alumnosSep[x] = students[x].sep;
    }

    ///////////////////////////////////////////////////
    /// Se asignan las escuelas un arreglo que y estudiantes a la escuela
    /// las escuelas tendran como identificación el indice
    /// y currentSolution tiene como indice al estudiante y el valor del indice a la escuela que asignada
    ///////////////////////////////////////////////////
    assignSchoolToArray(previousSolution, bestSolution, currentSolution, ptr_colegios, ptr_students, cupoArray);
    calcDist(ptr_colegios, ptr_students, distMat);

    ///////////////////////////////////////////////////
    /// Termina La fase de recolección de datos.
    /// Es necesario crear una funcion que empareje al estudiante con la escuela correspondiente segun su puesto
    /// en el arreglo ejemplo el 5 estudiante tiene rbd 4566 ese apunta al colegio que esta en la posicion
    /// 20 entonces cambio el rbd del estudiante a 20
    ///////////////////////////////////////////////////

    ///////////////////////////////////////////////////
    ///  Configuración de parametros
    ///////////////////////////////////////////////////

    double  costBestSolution,
            costPreviousSolution,
            costCurrentSolution,
            *ptr_alpha = &alpha[0],
            sumaAlpha = 0;

    ///////////////////////////////////////////////////
    /// Calcula el valor de los alpha
    ///////////////////////////////////////////////////
    for(x=0; x<3; x++){
        sumaAlpha +=alpha[x];
    }
    for(x=0; x<3; x++){
        alpha[x]= alpha[x]/(double)sumaAlpha;
    }


    ////////////////////////////////////////////////
    ////// Hace una calculo de rango de los promedios de las distancias
    ///////////////////////////////////////////////////

    for(int i=0;i<n_students;i++){
        for(x=0;x<n_colegios;x++){
            if(distMat[i][x]>max_dist){
                max_dist = distMat[i][x];
            }
        }
    }

    std::cout << "Max promedio dist:" << max_dist << "| Min Promedio dist: " << min_dist << "| init_dist: "<< init_dist<<"\n";


    ///////////////////////////////////////////////////
    /// Registro de datos
    ///////////////////////////////////////////////////

    costBestSolution=calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);
    std::cout << "--------------- Primeros datos -------------" << "\n";
    info << "--------------- Primeros datos -------------" << "\n";
    std::cout << "Primer costo de solución: " << costBestSolution << "\n";
    info << "Primer costo de solución: " << costBestSolution << "\n";
    costPreviousSolution=costBestSolution;
    costCurrentSolution=costBestSolution;
    std::cout << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    info << "Primer distancia: " << meanDist(currentSolution,distMat) << "\n";
    std::cout << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    info << "Primer Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
    std::cout << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";
    info << "Primer CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";


    ///////////////////////////////////////////////////
    /// Generación de archivos que almacenan información de los graficos
    ///////////////////////////////////////////////////
    std::ofstream info_graficos;
    std::string name_info_graficos = "./save/" + std::string(timestr) +"-info-graficos.txt"; // concatenar
    info_graficos.open(name_info_graficos);
    info_graficos << count << "," << meanDist(currentSolution,distMat) << "," << S(currentSolution, alumnosSep, totalVuln) << "," << costCupo(currentSolution,cupoArray) << "," << costCurrentSolution << "," << temp << "\n";

    ///////////////////////////////////////////////////
    /// Genera arreglos que contendran valores del 0 hasta n_students y n_colegios
    ///////////////////////////////////////////////////
    int *shuffle_student = new int[n_students];
    int *shuffle_colegios = new int[n_colegios];
    for (int i = 0; i < n_students; i++) {
        shuffle_student[i] = i;
    }
    for (int i=0; i < n_colegios; i++){
        shuffle_colegios[i]=i;
    }
    ///////////////////////////////////////////////////
    /// Posicion estudiantes
    ///////////////////////////////////////////////////

    std::ofstream info_graficos_bestSolution;
    std::string name_info_graficos_bestSolution = ruta_save + prefijo_save +"-info-graficos_bestSolution.txt"; // concatenar
    info_graficos_bestSolution.open(name_info_graficos_bestSolution);
    for(x=0;x<n_students;x++){
        info_graficos_bestSolution << currentSolution[x] << ",";
    }
    info_graficos_bestSolution << "\n";


    ///////////////////////////////////////////////////
    /// Genera distribuciones para seleccionar un estudiante y una escuela al azar
    ///////////////////////////////////////////////////

    dist = std::uniform_int_distribution<int>(0, n_students-1);
    dist2 = std::uniform_int_distribution<int>(0, n_colegios-1);


    ///////////////////////////////////////////////////
    /// Genera los archivos que contienen información de los estados de estudiantes y escuelas durante
    /// la ejecución del algoritmo
    ///////////////////////////////////////////////////
    std::ofstream studentscsv,schoolcsv;
    std::string nameCsvStudent = "./save/" + std::string(timestr) +"-students.csv"; // concatenar
    studentscsv.open(nameCsvStudent);
    studentscsv << "school,sep,latitude,longitude,count\n";
    std::string nameCsvSchool = "./save/" + std::string(timestr) +"-school.csv"; // concatenar
    schoolcsv.open(nameCsvSchool);
    schoolcsv << "id,latitude,longitude,rbd\n";

    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////
    auto start = std::chrono::high_resolution_clock::now();
    int count_rechaso=0;

    int count_acepta_cost=0;
    int count_acepta_suerte=0;
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////


    std::vector<double> vector_costCurrentSolution;
    std::vector<double> vector_meanDist;
    std::vector<double> vector_segregation;
    std::vector<double> vector_costoCupo;
    std::vector<double> vector_temp;
    std::vector<int> vector_count;
    int c_accepta = 0;
    count++;

    while(temp > min_temp){


        // Obtiene el costo Actual
        newSolution(currentSolution,previousSolution);
        costCurrentSolution = calCosto(currentSolution,distMat,ptr_alpha, alumnosSep, totalVuln, cupoArray);


        if(costCurrentSolution<0.00){

            std::cout << "distancia: " << meanDist(currentSolution,distMat) << "\n";
            std::cout << "Segregación: " << S(currentSolution, alumnosSep, totalVuln) << "\n";
            std::cout << "CostoCupo: " << costCupo(currentSolution,cupoArray) << "\n";
            std::cout << costCurrentSolution;
            exit(1);
        }
        // Verifica si el costo actual es mejor que la mejor solución
        // en el caso que el costo actual es menor a la mejor solución, acepta la solución y los
        // guarda en el estado como mejor solución
        //std::cout << "CostoCurreent segundo" << costCurrentSolution << "\n";




        if(costCurrentSolution < costBestSolution){
            // guarda la actual solución como la mejor
            for(x=0;x<n_students;x++){
                bestSolution[x]=currentSolution[x];
                previousSolution[x]=currentSolution[x];
            }

            costBestSolution=costCurrentSolution;
            costPreviousSolution=costCurrentSolution;

            vector_costCurrentSolution.push_back(costCurrentSolution);
            vector_meanDist.push_back(meanDist(currentSolution,distMat));
            vector_segregation.push_back(S(currentSolution, alumnosSep, totalVuln));
            vector_costoCupo.push_back(costCupo(currentSolution,cupoArray));
            vector_temp.push_back(temp);
            vector_count.push_back(count);
            c_accepta++;
            count_rechaso=0;
        }
        // En el caso que el la solución actual sea mas alta intenta aceptar una peor solución en base
        // a la función acepta
        else{
            if(acepta(costPreviousSolution,costCurrentSolution)==1){
                for(x=0;x<n_students;x++){
                    previousSolution[x]=currentSolution[x];
                }
                costPreviousSolution=costCurrentSolution;
                count_rechaso=0;
                c_accepta++;
            }
            else{
                count_rechaso++;
            }
        }


        ///////////////////////////////////////////////////
        /// Largo de temperatura
        ///////////////////////////////////////////////////
        if(c_accepta>=n_colegios){
            temp=temp*(coolingRate);
            //std::cout << "Enfriamiento " << temp << " CostZ " << costCurrentSolution << " bestZ " << costBestSolution << " count_rechaso " << count_rechaso << " count_reheating " << count_reheating <<"\n";
            c_accepta=0;
        }
        if(count%((n_colegios*2))==0){
            temp=temp*(coolingRate);
            //std::cout << "recalentamiento " << temp << " CostZ " << costCurrentSolution << " bestZ " << costBestSolution << " count_rechaso " << count_rechaso << " count_reheating " << count_reheating <<"\n";
        }

        ///////////////////////////////////////////////////
        /// Reinicio de temperatura
        ///////////////////////////////////////////////////
        count++;

    }

    std::cout << "cuenta acepta costo: " << count_acepta_cost << "\n";
    std::cout << "cuenta acepta suerte: " << count_acepta_suerte << "\n";

    ///////////////////////////////////////////////////
    /// Obtiene el tiempo de ejecución
    ///////////////////////////////////////////////////
    auto end = std::chrono::high_resolution_clock::now();
    double time_taken = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    time_taken *= 1e-9;

    for(x=0;x<n_students;x++){
        info_graficos_bestSolution << bestSolution[x] << ",";
    }
    info_graficos_bestSolution.close();

    for(x=0; x<vector_count.size(); x++){
        info_graficos << vector_count.at(x) << "," << vector_meanDist.at(x) << "," << vector_segregation.at(x) << "," << vector_costoCupo.at(x) << "," << vector_costCurrentSolution.at(x) << "," << std::fixed << vector_temp.at(x) << std::setprecision(13) << "\n";
    }
    info_graficos.close();



    ///////////////////////////////////////////////////
    /// Almacenamiento de datos
    ///////////////////////////////////////////////////

    std::cout << "--------------- Resultado Final ----------------" << "\n";
    info  << "--------------- Resultado Final ----------------" << "\n";
    std::cout << "Numero de Ciclos " << count << "\n";
    info  << "Numero de Ciclos " << count << "\n";
    std::cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    info << "Costo de la solución previa: " << costPreviousSolution << "\n";
    std::cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    info << "Costo de la mejor solución: " << costBestSolution << "\n";
    std::cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    info << "Costo de la solución actual: " << costCurrentSolution << "\n";
    std::cout << "Tiempo de ejecución de SA: " << std::fixed << time_taken << std::setprecision(9) << "\n";
    info << "Tiempo de ejecución de SA: " << std::fixed << time_taken << std::setprecision(9) << "\n";
    std::cout << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    info << "distancia: " << meanDist(bestSolution,distMat) << "\n";
    std::cout << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    info << "Segregación: " << S(bestSolution, alumnosSep, totalVuln) << "\n";
    std::cout << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";
    info << "CostoCupo: " << costCupo(bestSolution,cupoArray) << "\n";
    std::cout << "--------------- Finalizo con exito ----------------" << "\n";
    info << "--------------- Finalizo con exito ----------------" << "\n";

    std::ofstream info_test;
    std::string nameinfo_test = ruta_save + prefijo_save+"-info_test.txt"; // concatenar
    info_test.open(nameinfo_test);
    info_test << std::fixed << time_taken << std::setprecision(9) << "," << costBestSolution << "," << meanDist(bestSolution,distMat) << "," << S(bestSolution, alumnosSep, totalVuln) << "," << costCupo(bestSolution,cupoArray) << "," << count << "," << std::fixed << temp << std::setprecision(13) << "," << min_temp << "," << coolingRate << "," << alpha1 << "," << alpha2 << "," << alpha3 << "," << seed << "\n";
    info_test.close();

    std::cout << "-------------- Guardando Archivos /cmake-build-dbug-save -----------------" << "\n";
    info.close();
    studentscsv.close();
    schoolcsv.close();
    std::cout << "-------------- Archivos Guardado ------------------" << "\n";
    return (EXIT_SUCCESS);

}

///////////////////////////////////////////////////
/// Calcula el costo
///////////////////////////////////////////////////
double calCosto(int currentSolution[], double **distMat, const double ptr_alpha[], int alumnosSep[], int totalVuln, int cupoArray[]){
    double var1 = meanDist(currentSolution,distMat)/max_dist;
    //std::cout << "distancia: " << var1 << "\n";
    double var2 = S(currentSolution, alumnosSep, totalVuln);
    //std::cout << "Segregación: " << var2 << "\n";
    double var3 = costCupo(currentSolution,cupoArray);
    //std::cout << "CostoCupo: " << var3 << "\n";
    return (double)((ptr_alpha[0]*var1)+(ptr_alpha[1]*var2)+(ptr_alpha[2]*var3));
}

///////////////////////////////////////////////////
/// Distancia promedio que recorreran los estudiantes
///////////////////////////////////////////////////
double meanDist(const int currentSolution[], double  **distMat){
    double sumDist=0;
    for(int i=0;i<n_students;i++){
        sumDist+=distMat[i][currentSolution[i]]; // distMat[estudiante][escuela]
    }
    double mean=sumDist/double(n_students);
    //std::cout << "Numero de estudiantes: " << n_student << "  |  Suma de distancias:" << sumDist << "\n";
    return mean;
}

///////////////////////////////////////////////////
/// Calcula segregación por duncan
///////////////////////////////////////////////////

double S(const int currentSolution[],const int alumnosSep[], int totalVuln){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++){
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++){
            if(currentSolution[a] == n){
                aluNoVulCol++;
                aluVulCol+=alumnosSep[a];
            }
        }
        if(aluNoVulCol>0){
            aluNoVulCol =aluNoVulCol - aluVulCol;
            totalSesc+=((double)1/2)*fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln)));
        }
    }
    return totalSesc;
}

///////////////////////////////////////////////////
/// Calcula el costo de tener los estudiantes en las escuelas
///////////////////////////////////////////////////

double costCupo(const int currentSolution[],const int cupoArray[]){
    double totalcostCupo = 0;
    int totalAluCol = 0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
            if(currentSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+=totalAluCol*fabs((cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
    }
    return (totalcostCupo/n_colegios);
}

///////////////////////////////////////////////////
/// Genera una nueva solución en donde asigna a un estudiante a una escuela
/// aleatoriamente
///////////////////////////////////////////////////

void newSolution(int currentSolution[],const int previousSolution[]){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, n_students);
    std::random_device rd2;
    std::mt19937 mt2(rd2());
    std::uniform_int_distribution<int> dist2(0, n_colegios);
    int selectStudent=dist(mt);
    int selectSchool = dist2(mt2);
    for(int x=0; x<n_students; x++){
        if(x == selectStudent) {
            currentSolution[x] = selectSchool;
        }
        else {
            currentSolution[x] = previousSolution[x];
        }
    }

}

///////////////////////////////////////////////////
/// Función de aceptación en base a mayor temperatura mayor probabilidad que acepte a una solución peor
/// en caso de menor temperatura menor probabibilidad que acepte una solución peor.
///////////////////////////////////////////////////

int acepta(double costPrevious, double costCurrent){;
    std::uniform_real_distribution<double> dist_accepta(0.0, 1.0);
    if(costCurrent < costPrevious){
        return 1;
    }
    else{
        double valor=p(costPrevious,costCurrent);
        double nrandom=dist_accepta(mt);
        if(nrandom<valor){
            return 1;
        }
        else{
            return 0;
        }
    }
}
double p(double costPrevious,double costCurrent){
    double po;
    po = exp(-(costCurrent-costPrevious)/((double)temp));
    return po;
}

///////////////////////////////////////////////////
/// Asigna a las soluciones la escuela actual Solo se utiliza al inicio
///////////////////////////////////////////////////
void assignSchoolToArray(int previousSolution[], int bestSolution[], int currentSolution[], Info_colegio *ptr_colegios, Info_alu *ptr_students, int cupoArray[]){
    Info_alu *ptr_aux = ptr_students;
    for(int x=0;x < n_colegios;x++){
        for(int y=0; y < n_students; y++){
            if(ptr_colegios->rbd == ptr_students->rbd){
                previousSolution[y] = x;
                bestSolution[y] = x;
                currentSolution[y] = x;
            }
            ptr_students++;

        }
        /*
         * cupoArray sera un arreglo que por indice es la escuela y su valor sera el cupo que posee esa escuela
         * se asume que las escuelas pueden tener sobre cupo.
         */

        cupoArray[x] = ptr_colegios->num_alu+ ((int)((ptr_colegios->num_alu*10)/100));

        ptr_students = ptr_aux;
        ptr_colegios++;
    }
}
///////////////////////////////////////////////////
/// Crea una matriz de distancia donde x es el estudiante, y es la escuela
///////////////////////////////////////////////////
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students, double **distMat){
    Info_colegio *ptr_aux = ptr_colegios;
    for(int x=0;x < n_students ;x++){
        for(int y=0; y < n_colegios; y++){
            distMat[x][y] = sqrt( pow((ptr_students->latitude - ptr_colegios->latitude),2)+pow((ptr_students->longitude - ptr_colegios->longitude),2))/1000;
            ptr_colegios++;

        }
        ptr_colegios = ptr_aux;
        ptr_students++;
    }
}
