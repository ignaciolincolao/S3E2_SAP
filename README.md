## Requerimientos

Las versiones sequential y multihilo se conectan directamente a mongodb por lo que es necesario tener una base de datos y montarlo desde ahi.

### Librerias requeridas

#### C++

-libbson
-libmongoc
-mongocxx
-bsoncxx

#### Python

-json
-utm

## Instrucciones

Primero es necesario generar los estudiantes por lo que se debe ejecutar lo siguiente:

    cd generate_alu_python
    chmod +x alus_generate.py
    ./alus_generate.py

    
Esto genera los estudiantes y escuelas en distintos formatos txt y json, en el formato txt tienen las siguiente estructura los archivos: 

| alumnos_lat_lon.txt | RBD | Latitud | Longitud | Sep |
| ------------ | --- | ------- | -------- | --- |
|              | 5668 | -38.71250725 | -72.65644069 | 1 |
    
| alumnos_utm.txt | RBD | UTM Easting | UTM Northing | Sep |
| ------------ | --- | ------- | -------- | --- |
|              | 5668 | 703764.0664963166 | 5712518.935511605 | 1 |
    
| colegios_lat_lon.txt | RBD | Latitud | Longitud | Alumnos | Prioritarios |
| ------------ | --- | ------- | -------- | --- | --- |
|              | 5668 | -38.28600077 | -72.56208477 | 275 | 266 |
    
| colegios_utm.txt | RBD | UTM Easting | UTM Northing | Alumnos | Prioritarios |
| ------------ | --- | ------- | -------- | --- | --- |
|              | 5668 | -38.71250725 | -72.65644069 | 275 | 266 |

Los archivos json tienen poseen las misma estructura, para importar a la base de datos pueden ejecutar el archivo load_databse.sh:

```
chmod +x load_database.sh
./load_database.sh
```

Con los datos generados y cargados a la base de datos funcionan sin problemas las siguientes versiones:

-sequential
-sequential_v2
-multihilo

Para las versiones de paralelizado es necesario colocar los archivos.txt al lado del ejecutable, que se encuentra en las carpetas cmake-build-debug:

```
cp ./generate_alu_python/alumnos_lat_lon.txt /<paralelizado_..>/cmake-build-debug/alumnos_lat_lon.txt
cp ./generate_alu_python/alumnos_lat_lon.txt /<paralelizado_..>/cmake-build-debug/colegios_lat_lon.txt
cp ./generate_alu_python/alumnos_lat_lon.txt /<paralelizado_..>/cmake-build-debug/alumnos_utm.txt
cp ./generate_alu_python/alumnos_lat_lon.txt /<paralelizado_..>/cmake-build-debug/colegios_utm.txt
```

    
Para compilar se debe entrar al cmake-build-debud de la version que se quiera ejecutar, una ves dentro ejecutar los siguientes comandos:

```
cmake ..
cmake --build .
```
## Información de los programas

- recocido_simulado_paralelizado: escoge el mejor resultado, vector de escuelas es el mismo en cada estudiante
- recocido_simulado_paralelizado_v2: escoge el primer mejor resultado vector de escuelas es el mismo en cada estudiante
- recocido simulado paralelizado_v3: escoge el mejor resultado, vector de escuelas es diferente para cada estudiante
- recocido simulado paralelizado_v4: escoge el mejor resultado, vector de escuelas es el mismo en cada estudiante (selección paralela corregida)


