#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <cuda.h>
#include <chrono>
#include "UTM.h"
#include <iomanip>

///////////////////////////////////////////////////
/// Variables constantes CUDA
///////////////////////////////////////////////////

__constant__ int d_cupoArray[85];
__constant__ double d_alpha[3];



///////////////////////////////////////////////////
/// Estructura de datos de los colegios.
///////////////////////////////////////////////////
struct Info_colegio {
    double latitude = 0;
    double longitude = 0;
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

double calCosto(const int currentThreadSolution[]);
double meanDist(const int currentThreadSolution[]);
double S(const int currentThreadSolution[]);
double costCupo(const int currentThreadSolution[]);
int acepta(double costPrevious, double costCurrent);
double p(double costPrevious,double costCurrent);
void assignSchoolToArray(Info_colegio *ptr_colegios, Info_alu *ptr_students);
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students);
void shuffle(int[],int,std::uniform_int_distribution<int>);

///////////////////////////////////////////////////
/// Kernel newSolution_kerneln_colegios
///////////////////////////////////////////////////


__global__ void newSolution_kernel(
        double *d_array_current_Solution,
        int *d_array_current_Solution_thread,
        const int n_students,
        const int n_colegios,
        const int n_thread,
        const double max_dist,
        const int* __restrict__ d_alumnosSep,
        const int totalVuln,
        const int* __restrict__ d_aluxcol,
        const int* __restrict__ d_aluVulxCol,
        const int* __restrict__ d_currentSolution,
        const double* __restrict__ d_distMat,
        const int* __restrict__ d_shuffle_students,
        const int* __restrict__ d_shuffle_colegios,
        size_t pitch){

    /// Shared Memory
    extern __shared__ double sharedMem[];
    int* aluxcolblock = (int*)sharedMem;
    int* aluVulxColblock = (int*)&aluxcolblock[n_colegios];
    double* solutions =(double*)&aluVulxColblock[n_colegios];
    int* solutions_thread = (int*)&solutions[n_thread];
    /// Inicializa variables en 0
    int aluchange,
            colchange,
            i = 0,
            x = 0,
            aluVulCol= 0,
            aluNoVulCol= 0,
            totalAluCol= 0,
            myID = threadIdx.x,
            school_alu_change,
            salto= n_thread;

    double  totalcostCupo= 0.0,
            totalSesc= 0.0,
            var1,
            var2,
            var3,
            result= 0.0;
    /// Inicializa arrays
    aluchange = d_shuffle_students[blockIdx.x];
    colchange = d_shuffle_colegios[threadIdx.x];
    solutions_thread[threadIdx.x] = colchange;

    /// Recopila la informacion que existe en memoria global
    clock_t start_time = clock();
    /// a shared memory29853
    school_alu_change = d_currentSolution[aluchange];
    for (i = threadIdx.x; i< n_colegios; i=i+n_thread){
        aluxcolblock[i] = d_aluxcol[i];
        aluVulxColblock[i] = d_aluVulxCol[i];
        if(i == school_alu_change){
            aluxcolblock[school_alu_change]-=1;
            aluVulxColblock[school_alu_change]-=d_alumnosSep[aluchange];
        }
    }

    /// Calcula la distancia total
    for (x = 0 ; x < n_students ; x++) {
        if (x != aluchange) {
            result += d_distMat[x * pitch / sizeof(double) + d_currentSolution[x]];
        }
        else {
            result += d_distMat[x * pitch / sizeof(double) + colchange];
        }
    }
    __syncthreads();
    clock_t stop_time = clock();
    int valtime = (int)(stop_time - start_time);
    if(threadIdx.x==0 && blockIdx.x==0){
        printf("%d \n",valtime);
    }
    /// Calcula el costo cupo y la cantidad de segregación total
    for(int n=0; n<n_colegios; n++){
        totalAluCol = aluxcolblock[n];
        aluVulCol = aluVulxColblock[n];
        if(n == colchange){
            totalAluCol+=1;
            aluVulCol+=d_alumnosSep[aluchange];
        }
        aluNoVulCol =totalAluCol - aluVulCol;
        // Calcula el costo cupo
        totalcostCupo+=totalAluCol*fabs((d_cupoArray[n]-totalAluCol)/pow(((double)d_cupoArray[n]/2),2));
        // Calcula el total sesc
        totalSesc+=((double)1/2)*fabs((aluVulCol/(double)totalVuln)-(aluNoVulCol/(double)(n_students-totalVuln)));
    }

    var1 = d_alpha[0]*((result/(double(n_students)))/double(max_dist));
    var2 = d_alpha[1]*totalSesc;
    var3 = d_alpha[2]*(totalcostCupo/n_colegios);
    solutions[myID] = var1+var2+var3;
    if(colchange == school_alu_change){
        solutions[myID] = 1.0;
    }
    __syncthreads();
    while(salto){
        if(salto-(myID+1)>myID){
            if(solutions[myID]>solutions[salto-(myID+1)]){
                solutions[myID]=solutions[salto-(myID+1)];
                solutions_thread[myID]=solutions_thread[salto-(myID+1)];
            }
        }
        salto = (salto/2)+(salto&(2-1));
        if(salto==1){
            salto = 0;
        }
        __syncthreads();
    }
    if(myID==0)
    {
        d_array_current_Solution[blockIdx.x] = solutions[myID];
        d_array_current_Solution_thread[blockIdx.x] = solutions_thread[myID];

    }
}

__global__ void reduce_block_kernel(
        double *d_array_current_Solution,
        int *d_array_current_Solution_thread,
        int *d_array_current_Solution_block,
        const int n_block){

    extern __shared__ double sharedMem[];
    double* solutions =(double*)sharedMem;
    int* solutions_block = (int*)&solutions[n_block];
    int* solutions_thread = (int*)&solutions_block[n_block];

    int myID = threadIdx.x;
    int salto= n_block;
    solutions[myID] = d_array_current_Solution[myID];
    solutions_thread[myID] = d_array_current_Solution_thread[myID];
    solutions_block[myID]= myID;
    __syncthreads();
    while(salto){
        if(salto-(myID+1)>myID){
            if(solutions[myID]>solutions[salto-(myID+1)]){
                solutions[myID]=solutions[salto-(myID+1)];
                solutions_thread[myID]=solutions_thread[salto-(myID+1)];
                solutions_block[myID]=solutions_block[salto-(myID+1)];
            }
        }
        salto = (salto/2)+(salto&(2-1));
        if(salto==1){
            salto = 0;
        }
        __syncthreads();
    }
    if(myID==0)
    {
        
        d_array_current_Solution[myID] = solutions[myID];
        d_array_current_Solution_thread[myID]= solutions_thread[myID];
        d_array_current_Solution_block[myID] = solutions_block[myID];
    }
}



///////////////////////////////////////////////////
/// Parametros de configuración Default
///////////////////////////////////////////////////

double alpha1 = 15; // Alpha de distancia
double alpha2 = 30; // Alpha de segregación
double alpha3 = 25; // Alpha de costocupo
double coolingRate = 0.98; // Tasa de enfriamiento
double temp = 100000; // Temperatura inicial
double min_temp = 0.00000009; // Minima temperatura que puede llegar
int n_block = 256; // Numero de blockes = numeros de alumnos aleatorios
int n_thread = 85; // Numero de threads por bloque = numeros de escuelas aleatorios
std::string ruta_save = "./save/"; // Ruta para guardar los archivos
double k_recalentamiento = 0.90;
double max_temp = 0;
double e_const=0.01;
int count_rechaso=0;


///////////////////////////////////////////////////
/// Variables globales.
///////////////////////////////////////////////////
double alpha[3]={alpha1,alpha2,alpha3},
        **distMat=nullptr;
int n_students = 0,
        n_colegios,
        selectThread=0,
        selectBlock = 0,
        totalVuln = 0,
        *alumnosSep=nullptr,
        *cupoArray=nullptr;
int *previousSolution= nullptr;
int *bestSolution= nullptr;
int *currentSolution=nullptr;
int seed= 280;//rand();
double max_dist=0.0;
std::random_device rd;
std::mt19937 mt(rd());
std::uniform_int_distribution<int> dist(0,0);
std::uniform_int_distribution<int> dist2(0,0);


///////////////////////////////////////////////////
/// Funcion principal
///////////////////////////////////////////////////

int main(int argc, char *argv[]) {
    int  test_vectosize[29853];
    int  test_vectosize2[85*3];
    std::cout << "peso del vector" << sizeof(test_vectosize2) << "\n";

    time_t hora_actual;
    struct tm *time_info;
    time(&hora_actual);
    time_info = localtime(&hora_actual);
    char timestr[20];
    strftime(timestr, sizeof(timestr), "%Y-%m-%d T:%H-%M", time_info);
    std::string prefijo_save = std::string(timestr);

    if (argc>1) {
        alpha1 = std::stod(argv[1]); // Alpha de distancia
        alpha2 = std::stod(argv[2]); // Alpha de segregación
        alpha3 = std::stod(argv[3]); // Alpha de costocupo
        alpha[0]=alpha1;
        alpha[1]=alpha2;
        alpha[2]=alpha3;
        coolingRate = std::stod(argv[4]); // Tasa de enfriamiento
        k_recalentamiento = std::stod(argv[5]);
        temp = std::stod(argv[6]); // Temperatura inicial
        min_temp = std::stod(argv[7]); // Minima temperatura que puede llegar
        n_block = std::stoi(argv[8]); // Numero de blockes = numeros de alumnos aleatorios
        n_thread = std::stoi(argv[9]); // Numero de threads por bloque = numeros de escuelas aleatorios
        ruta_save = argv[10];
        prefijo_save = argv[11];
        max_temp= pow(10,300);
        seed= std::stoi(argv[12]);
    }
    mt.seed(seed);
    Info_colegio *ptr_colegios;
    Info_alu *ptr_students;




    std::ofstream info;
    std::string infotxt = ruta_save + prefijo_save + "-info.txt"; // concatenar
    info.open(infotxt);
    int x = 0, z = 0;

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
    /// Se crea las soluciones que tienen n_students de largo
    ///////////////////////////////////////////////////
    previousSolution = (int *) malloc(sizeof(int) * n_students);
    bestSolution = (int *) malloc(sizeof(int) * n_students);



    ///////////////////////////////////////////////////
    /// Se asignan las escuelas un arreglo que y estudiantes a la escuela
    /// las escuelas tendran como identificación el indice
    /// y currentSolution tiene como indice al estudiante y el valor del indice a la escuela que asignada
    ///////////////////////////////////////////////////

    currentSolution = (int *) malloc(sizeof(int) * n_students);
    /// Se crea una matriz de distnacia donde se obtienen todas las distancias entre estudiantes y escuelas.
    distMat = (double **) malloc(sizeof(double) * n_students);
    cupoArray = (int *) malloc(sizeof(int) * n_colegios);
    ///Alumnos sep
    alumnosSep = (int *) malloc(sizeof(int) * n_students);
    /// Se crear un arreglo donde el el valor es la posición del estudiante sep
    for (x = 0; x < n_students; x++) {
        distMat[x] = (double *) malloc(sizeof(double) * n_colegios);
        alumnosSep[x] = students[x].sep;
    }

    assignSchoolToArray(ptr_colegios, ptr_students);
    calcDist(ptr_colegios, ptr_students);


    ///////////////////////////////////////////////////
    /// Termina La fase de recolección de datos.
    /// Es necesario crear una funcion que empareje al estudiante con la escuela correspondiente segun su puesto
    /// en el arreglo ejemplo el 5 estudiante tiene rbd 4566 ese apunta al colegio que esta en la posicion
    /// 20 entonces cambio el rbd del estudiante a 20
    ///////////////////////////////////////////////////

    double costBestSolution,
            costPreviousSolution,
            costCurrentSolution,
            sumaAlpha = 0;

    ///////////////////////////////////////////////////
    /// Calcula el valor de los alpha
    ///////////////////////////////////////////////////

    for (x = 0; x < 3; x++) {
        sumaAlpha += alpha[x];
    }

    for (x = 0; x < 3; x++) {
        alpha[x] = alpha[x] / (double) sumaAlpha;
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


    ///////////////////////////////////////////////////
    /// Registro de datos
    ///////////////////////////////////////////////////
    costBestSolution = calCosto(currentSolution);
    std::cout << "Primer costo de solución: " << costBestSolution << "\n";
    costPreviousSolution = costBestSolution;
    costCurrentSolution = costBestSolution;
    std::cout << "Primer distancia: " << meanDist(currentSolution) << "\n";


    std::cout << "Primer Segregación: " << S(currentSolution) << "\n";
    std::cout << "Primer CostoCupo: " << costCupo(currentSolution) << "\n";
    int count = 0;

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    /// Reservando memoria en CUDA
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    ///////////////////////////////////////////////////
    /// Matrices 2D
    ///////////////////////////////////////////////////
    double *d_distMat; /// clon de matriz de distancia
    int *d_currentSolution;
    int *d_alumnosSep; // Array que contendra a los estudiantes vulnerables

    ///////////////////////////////////////////////////
    /// Valores para las partes calcular la nueva solucion
    ///////////////////////////////////////////////////

    int aluVulxCol[n_colegios], aluxcol[n_colegios];
    int previousAluxCol[n_colegios];
    int previousAluVulxCol[n_colegios];
    int bestAluxCol[n_colegios];
    int bestAluVulxCol[n_colegios];

    for(x = 0; x < n_colegios; x++){
        aluxcol[x] = colegios[x].num_alu;
        previousAluxCol[x] = colegios[x].num_alu;
        bestAluxCol[x] = colegios[x].num_alu;
        aluVulxCol[x] = colegios[x].prioritario;
        previousAluVulxCol[x] = colegios[x].prioritario;
        bestAluVulxCol[x] = colegios[x].prioritario;

    }

    ///////////////
    double *d_array_current_Solution;
    int *d_array_current_Solution_thread;
    int *d_array_current_Solution_block;
    ///////////////
    int *d_aluxcol;
    int *d_aluVulxCol;
    int *d_shuffle_students;
    int *d_shuffle_colegios;




    cudaMalloc((void **) &d_array_current_Solution, n_block * sizeof(double));
    cudaMalloc((void **) &d_array_current_Solution_thread, n_block * sizeof(int));
    cudaMalloc((void **) &d_array_current_Solution_block, sizeof(int));
    cudaMalloc((void **) &d_shuffle_colegios, n_thread * sizeof(int));
    cudaMalloc((void **) &d_shuffle_students, n_block * sizeof(int));
    cudaMalloc((void **) &d_aluxcol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_aluVulxCol,n_colegios * sizeof(int));
    cudaMalloc((void **) &d_currentSolution, n_students * sizeof(int));  // Solución actual
    cudaMalloc((void **) &d_alumnosSep, n_students * sizeof(int)); // arreglo que contiene la id de cada usuario vulnerable
    double *matrestest = (double *) malloc(sizeof(double) * n_students * n_colegios);
    double *array_costCurrentSolution = (double *) malloc(sizeof(double) * n_block * n_thread);
    for (x = 0; x < n_students; x++) {
        for (z = 0; z < n_colegios; z++) {
            matrestest[n_colegios * x + z] = distMat[x][z];
        }
    }
    for (x = 0; x < n_block; x++){
        for (z = 0; z < n_thread; z++){
            array_costCurrentSolution[n_thread * x + z] = 0.0;
        }
    }





    ///////////////////////////////////////////////////
    /// Valores que nunca van a cambiar
    //////////////////////////////////////////////////////

    cudaMemcpy(d_alumnosSep, alumnosSep, n_students * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol( d_cupoArray, cupoArray,  n_colegios * sizeof(int));
    cudaMemcpyToSymbol( d_alpha, alpha,  3 * sizeof(double));

    size_t pitch;
    cudaMallocPitch(&d_distMat,
                    &pitch,
                    n_colegios * sizeof(double),
                    n_students); // Reserva memoria para la matriz de distancia

    size_t h_pitchBytes = n_colegios * sizeof(double);
    cudaMemcpy2D(d_distMat,
                 pitch,
                 matrestest,
                 h_pitchBytes,
                 n_colegios * sizeof(double),
                 n_students,
                 cudaMemcpyHostToDevice);


    ///////////////////////////////////////////////////
    /// Genera los archivos para graficos
    ///////////////////////////////////////////////////

    std::ofstream info_graficos;
    std::string name_info_graficos = ruta_save + prefijo_save +"-info-graficos.txt"; // concatenar
    info_graficos.open(name_info_graficos);
    info_graficos << count << "," << meanDist(currentSolution) << "," << S(currentSolution) << "," << costCupo(currentSolution) << "," << costCurrentSolution << "," << std::fixed <<  temp << std::setprecision(13) << "\n";
    std::cout << "Numero de bloques: " << n_block << "| Numeros de thread: "<< n_thread <<  "\n";
    strftime(timestr, sizeof(timestr), "%T", time_info);
    info << "Tiempo antes iniciar el ciclo" << prefijo_save << "\n";

    /// Posicion estudiantes

    std::ofstream info_graficos_bestSolution;
    std::string name_info_graficos_bestSolution = ruta_save + prefijo_save +"-info-graficos_bestSolution.txt"; // concatenar
    info_graficos_bestSolution.open(name_info_graficos_bestSolution);
    for(x=0;x<n_students;x++){
        info_graficos_bestSolution << currentSolution[x] << ",";
    }
    info_graficos_bestSolution << "\n";



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
    /// Inicializa las distribuciónes
    ///////////////////////////////////////////////////

    dist = std::uniform_int_distribution<int>(0, n_students-1);
    dist2 = std::uniform_int_distribution<int>(0, n_colegios-1);

    ///////////////////////////////////////////////////
    /// Contador de tiempo de ejecución en cuda
    ///////////////////////////////////////////////////

    cudaEvent_t start_cuda;
    cudaEvent_t stop_cuda;
    cudaEventCreate(&start_cuda);
    cudaEventCreate(&stop_cuda);
    float elapsedTime;

    ///////////////////////////////////////////////////
    /// Inicio el contador de tiempo antes de iniciar el algortimo
    ///////////////////////////////////////////////////


    cudaMemcpy(d_currentSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_aluxcol, aluxcol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_aluVulxCol, aluVulxCol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);





    std::vector<double> vector_costCurrentSolution;
    std::vector<double> vector_meanDist;
    std::vector<double> vector_segregation;
    std::vector<double> vector_costoCupo;
    std::vector<double> vector_temp;
    std::vector<int> vector_count;
    int reheating = 0;
    int c_accepta = 0;
    count++;
    double timeCuda = 0.0;

    int valmaxheating=25;
    int count_reheating = 0;
    double bestTemp = 0;

    double copytimecuda = 0;
    double executiontimecuda = 0;
    double recoverdatacuda = 0;
    ///////////////////////////////////////////////////
    /// Comienza a ejecutarse el algoritmo de SA
    ///////////////////////////////////////////////////
    auto start = std::chrono::high_resolution_clock::now();

    while(temp > min_temp){

        for(x=0;x<n_students;x++){
            currentSolution[x]=previousSolution[x];
        }
        for(x = 0; x < n_colegios; x++){
            aluxcol[x]=previousAluxCol[x];
            aluVulxCol[x]=previousAluVulxCol[x];
        }

        ///////////////////////////////////////////////////
        ///  Selecciona aleatoria mente a los alumnos
        ///////////////////////////////////////////////////

        shuffle(shuffle_student,n_block,dist);
        shuffle(shuffle_colegios,n_thread,dist2);

        ///////////////////////////////////////////////////
        /// Actualiza la memoria en CUDA
        ///////////////////////////////////////////////////


        cudaMemcpy(d_currentSolution, currentSolution, n_students * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_aluxcol, aluxcol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_aluVulxCol, aluVulxCol, n_colegios * sizeof(int), cudaMemcpyHostToDevice);




        ///////////////////////////////////////////////////
        ///  Envia datos a GPU
        ///////////////////////////////////////////////////

        cudaEventRecord(start_cuda,0);
        cudaMemcpy(d_shuffle_students, shuffle_student, n_block * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_shuffle_colegios, shuffle_colegios, n_thread * sizeof(int), cudaMemcpyHostToDevice);

        cudaDeviceSynchronize();
        cudaEventRecord(stop_cuda,0);
        cudaEventSynchronize(stop_cuda);
        cudaEventElapsedTime(&elapsedTime,start_cuda,stop_cuda);
        copytimecuda = copytimecuda+elapsedTime;

        ///////////////////////////////////////////////////
        ///  Ejecuta los kernel


        cudaEventRecord(start_cuda,0);
        newSolution_kernel<<<n_block,n_thread,
                n_colegios* sizeof(int) + n_colegios* sizeof(int) + n_thread* sizeof(double)+ n_thread* sizeof(int)>>>(
                        d_array_current_Solution,
                                d_array_current_Solution_thread,
                                n_students,
                                n_colegios,
                                n_thread,
                                max_dist,
                                d_alumnosSep,
                                totalVuln,
                                d_aluxcol,
                                d_aluVulxCol,
                                d_currentSolution,
                                d_distMat,
                                d_shuffle_students,
                                d_shuffle_colegios,
                                pitch);
        cudaDeviceSynchronize();

        reduce_block_kernel<<<1,n_block,
                n_block* sizeof(double)+ n_block* sizeof(int)+ n_block* sizeof(int)>>>(d_array_current_Solution,
                        d_array_current_Solution_thread,
                        d_array_current_Solution_block,
                        n_block);
        cudaDeviceSynchronize();

        cudaEventRecord(stop_cuda,0);
        cudaEventSynchronize(stop_cuda);
        cudaEventElapsedTime(&elapsedTime,start_cuda,stop_cuda);
        executiontimecuda = executiontimecuda+elapsedTime;

        cudaDeviceSynchronize();
        cudaEventRecord(start_cuda,0);
        cudaMemcpy(&costCurrentSolution,&d_array_current_Solution[0], sizeof(double),cudaMemcpyDeviceToHost);
        cudaMemcpy(&selectThread,&d_array_current_Solution_thread[0], sizeof(int),cudaMemcpyDeviceToHost);
        cudaMemcpy(&selectBlock,d_array_current_Solution_block, sizeof(int),cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        cudaEventRecord(stop_cuda,0);
        cudaEventSynchronize(stop_cuda);
        cudaEventElapsedTime(&elapsedTime,start_cuda,stop_cuda);
        recoverdatacuda = recoverdatacuda+elapsedTime;

        ///////////////////////////////////////////////////
        ///  Actualizo datos basicos
        ///////////////////////////////////////////////////

        //std::cout << "CPU: " << costCurrentSolution << " | " << selectThread << " | " << selectBlock << "\n";
        aluxcol[currentSolution[shuffle_student[selectBlock]]]-=1; ///
        aluVulxCol[currentSolution[shuffle_student[selectBlock]]]-=alumnosSep[shuffle_student[selectBlock]]; ///

        aluxcol[selectThread]+=1; ///
        aluVulxCol[selectThread]+=alumnosSep[shuffle_student[selectBlock]]; ///
        currentSolution[shuffle_student[selectBlock]] = selectThread; ///
        //aluxcol[shuffle_colegios[selectThread]]+=1; ///
        //aluVulxCol[shuffle_colegios[selectThread]]+=alumnosSep[shuffle_student[selectBlock]]; ///
        //currentSolution[shuffle_student[selectBlock]] = shuffle_colegios[selectThread]; ///




        ///////////////////////////////////////////////////
        /// Salida en caso de error
        ///////////////////////////////////////////////////
        //std::cout << costCurrentSolution << "\n";

        if(costCurrentSolution<0.00){
            std::cout << shuffle_colegios[selectThread] << "\n";
            std::cout << shuffle_student[selectBlock] << "\n";
            std::cout << "distancia: " << meanDist(currentSolution) << "\n";
            std::cout << "Segregación: " << S(currentSolution) << "\n";
            std::cout << "CostoCupo: " << costCupo(currentSolution) << "\n";
            std::cout << costCurrentSolution;
            exit(1);
        }






        if(costCurrentSolution < costBestSolution){
            // guarda la actual solución como la mejor
            for(x=0;x<n_students;x++){
                bestSolution[x]=currentSolution[x];
                previousSolution[x]=currentSolution[x];
            }
            for(x = 0; x < n_colegios; x++){
                previousAluxCol[x] = aluxcol[x];
                bestAluxCol[x] = aluxcol[x];
                previousAluVulxCol[x] = aluVulxCol[x];
                bestAluVulxCol[x] = aluVulxCol[x];

            }
            costBestSolution=costCurrentSolution;
            costPreviousSolution=costCurrentSolution;



            vector_costCurrentSolution.push_back(costCurrentSolution);
            vector_meanDist.push_back(meanDist(currentSolution));
            vector_segregation.push_back(S(currentSolution));
            vector_costoCupo.push_back(costCupo(currentSolution));
            vector_temp.push_back(temp);
            vector_count.push_back(count);
            c_accepta++;
            count_rechaso=0;
        }
        else{
            if(acepta(costPreviousSolution,costCurrentSolution)==1){
                /// Si por al asar acepta tomara la solución actual como la nueva solución a seguir
                for(x=0;x<n_students;x++){
                    previousSolution[x]=currentSolution[x];
                }
                for(x = 0; x < n_colegios; x++){
                    previousAluxCol[x] = aluxcol[x];
                    previousAluVulxCol[x] = aluVulxCol[x];
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
            //std::cout << "Enfriamiento " << temp << " CostZ " << costCurrentSolution  << " bestZ " <<  costBestSolution << " count_rechaso " << count_rechaso << " count_reheating " << count_reheating <<"\n";
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
        info_graficos << vector_count.at(x) << "," << vector_meanDist.at(x) << "," << vector_segregation.at(x) << "," << vector_costoCupo.at(x) << "," << vector_costCurrentSolution.at(x) << "," << std::fixed << vector_temp.at(x) << std::setprecision(13) << ","<< max_dist << "\n";
    }
    info_graficos.close();
    info << "Tiempo despues de finalizar el ciclo" << prefijo_save << "\n";
    std::cout << "--------------- Resultado Final ----------------" << "\n";
    info  << "--------------- Resultado Final con restauración----------------" << "\n";
    std::cout << "Numero de Ciclos: " << count << "\n";
    info << "Numero de Ciclos: " << count << "\n";
    std::cout << "Costo de la solución previa: " << costPreviousSolution << "\n";
    info << "Costo de la solución previa: " << costPreviousSolution << "\n";
    std::cout << "Costo de la mejor solución: " << costBestSolution << "\n";
    info << "Costo de la mejor solución: " << costBestSolution << "\n";
    std::cout << "Costo de la solución actual: " << costCurrentSolution << "\n";
    info << "Costo de la solución actual: " << costCurrentSolution << "\n";
    std::cout << "Tiempo de ejecución de SA: " << std::fixed << time_taken << std::setprecision(9) << "\n";
    info << "Tiempo de ejecución de SA: " << std::fixed << time_taken << std::setprecision(9) << "\n";
    std::cout << "distancia: " << meanDist(bestSolution) << "\n";
    info << "distancia: " << meanDist(bestSolution) << "\n";
    std::cout << "Segregación: " << S(bestSolution) << "\n";
    info << "Segregación: " << S(bestSolution) << "\n";
    std::cout << "CostoCupo: " << costCupo(bestSolution) << "\n";
    info << "CostoCupo: " << costCupo(bestSolution) << "\n";
    std::cout << "distancia * alpha: " << meanDist(bestSolution) * alpha[0] << "\n";
    info << "distancia * alpha: " << meanDist(bestSolution) * alpha[0] << "\n";
    std::cout << "Segregación * alpha: " << S(bestSolution) * alpha[1] << "\n";
    info << "Segregación * alpha: " << S(bestSolution) * alpha[1] << "\n";
    std::cout << "CostoCupo * alpha: " << costCupo(bestSolution) * alpha[2] << "\n";
    info << "CostoCupo * alpha: " << costCupo(bestSolution) * alpha[2] << "\n";
    std::cout << "solución calculada aparte: " << calCosto(bestSolution) << "\n";
    info << "solución calculada aparte: " << calCosto(bestSolution) << "\n";
    info << "Numero de bloques: " << n_block << "| Numeros de thread: "<< n_thread <<  "\n";
    std::cout << "Numero de bloques: " << n_block << "| Numeros de thread: "<< n_thread <<  "\n";
    std::cout << "--------------- Finalizo con exito ----------------" << "\n";
    info << "--------------- Finalizo con exito ----------------" << "\n";
    info.close();

    std::cout << "-------------- Guardando Archivos /cmake-build-dbug-save -----------------" << "\n";



    std::ofstream studentscsv,schoolcsv, info_test;

    std::string nameinfo_test = ruta_save + prefijo_save +"-info_test.txt"; // concatenar
    info_test.open(nameinfo_test);
    info_test << std::fixed << time_taken << std::setprecision(9) << "," << costBestSolution << "," << meanDist(bestSolution) << "," << S(bestSolution) << "," << costCupo(bestSolution) << ","<< n_block << "," << n_thread << "," << count << "," << std::fixed << temp << std::setprecision(13) << "," << min_temp << "," << coolingRate << "," << k_recalentamiento << "," << alpha1 << "," << alpha2 << "," << alpha3 << "," << seed <<  "," << timeCuda/1000 << "\n";
    info_test.close();

    std::cout << "Tiempo de copia: "<< copytimecuda << "\n";
    std::cout << "Tiempo de proceso: " <<  executiontimecuda << "\n";
    std::cout << "Tiempo de rescate de datos: " << recoverdatacuda << "\n";
    std::cout << "-------------- Archivos Guardado ------------------" << "\n";

    cudaFree(d_currentSolution);
    cudaFree(d_alumnosSep);
    cudaFree(d_distMat);
    cudaFree(d_cupoArray);
    cudaFree(d_array_current_Solution);
    cudaFree(d_array_current_Solution_thread);
    cudaFree(d_array_current_Solution_block);
    cudaFree(d_alpha);
    cudaEventDestroy(start_cuda);
    cudaEventDestroy(stop_cuda);
    return (EXIT_SUCCESS);

}

///////////////////////////////////////////////////
/// Calcula el costo
///////////////////////////////////////////////////
double calCosto(const int currentThreadSolution[]){
    double var1 = meanDist(currentThreadSolution)/max_dist;
    std::cout << "distancia: " << var1 << "\n";
    double var2 = S(currentThreadSolution);
    std::cout << "Segregación: " << var2 << "\n";
    double var3 = costCupo(currentThreadSolution);
    std::cout << "CostoCupo: " << var3 << "\n";
    return (double)((alpha[0]*var1)+(alpha[1]*var2)+(alpha[2]*var3));
}

///////////////////////////////////////////////////
/// Distancia promedio que recorreran los estudiantes
///////////////////////////////////////////////////
double meanDist(const int currentThreadSolution[]){
    double sumDist=0;
    for(int i=0;i<n_students;i++){
        sumDist+=distMat[i][currentThreadSolution[i]]; // distMat[estudiante][escuela]
    }
    double mean=sumDist/double(n_students);

    //std::cout << "Numero de estudiantes: " << n_student << "  |  Suma de distancias:" << sumDist << "\n";
    return mean;
}
///////////////////////////////////////////////////
/// Calcula segregación por duncan
///////////////////////////////////////////////////
double S(const int currentThreadSolution[]){
    double totalSesc = 0.0;
    int aluVulCol =0;
    int aluNoVulCol = 0;
    for(int n=0; n<n_colegios;n++) {
        aluVulCol = 0;
        aluNoVulCol = 0;
        for (int a = 0; a < n_students; a++) {
            if (currentThreadSolution[a] == n) {
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
double costCupo(const int currentThreadSolution[]){
    double totalcostCupo = 0;
    int totalAluCol = 0;
    for(int j=0;j<n_colegios;j++){
        totalAluCol = 0;
        for(int i=0; i<n_students; i++){
            if(currentThreadSolution[i]==j){
                totalAluCol++;
            }
        }
        totalcostCupo+=totalAluCol*fabs((cupoArray[j]-totalAluCol)/pow(((double)cupoArray[j]/2),2));
    }
    return (totalcostCupo/n_colegios);
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
    //po = 1/(1+exp(-(costPrevious-costCurrent)/temp));
    return po;
}

///////////////////////////////////////////////////
/// Asigna a las soluciones la escuela actual Solo se utiliza al inicio
///////////////////////////////////////////////////
void assignSchoolToArray(Info_colegio *ptr_colegios, Info_alu *ptr_students){
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
        //std::cout << ptr_colegios->num_alu+3 << "\n";
        cupoArray[x] = ptr_colegios->num_alu+ ((int)((ptr_colegios->num_alu*10)/100));
        ptr_students = ptr_aux;
        ptr_colegios++;
    }
}

///////////////////////////////////////////////////
/// Crea una matriz de distancia donde x es el estudiante, y es la escuela
///////////////////////////////////////////////////
void calcDist(Info_colegio *ptr_colegios, Info_alu *ptr_students){
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

void shuffle(int values[], const int max_change, std::uniform_int_distribution<int> distri) {
    int randvalue1,randvalue2,tem_value;
    for (int i = 0; i<max_change; i++) {
        randvalue1 = distri(mt);
        randvalue2 = i;
        tem_value = values[randvalue1];
        values[randvalue1] = values[randvalue2];
        values[randvalue2] = tem_value;
    }
    
}


///

//////// convierte la distancia
//// exp(x*fabs(((x-init)/max_dist)-init)/max_dist-min_dist)*exp(-(max_dist*fabs(((max_dist-init)/max_dist)-init)/max_dist-min_dist))