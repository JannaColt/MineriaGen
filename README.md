# MINERÍA GENÓMICA 
## CALIDAD, ENSAMBLE Y ANOTACIÓN
Contiene pipeline para ensamble y minería genómica de shotgun sequence con google colab

Desde el cuaderno establecido primeramente se instalan todos los paquetes que se usarán y al final se monta el drive en el que se estará trabajando. Es preferible que esto se haga desde el inicio ya que cuando se instala un nuevo paquete se reinicia el entorno y lo que anteriormente llamamos ya no estará disponible. 

## 1. Instalación de Herramientas

Instalamos conda y lo llamamos para proceder con la instalación de los demás paquetes usando conda
```python
!pip install -q condacolab
import condacolab
condacolab.install()
```
Luego instalamos miniconda, para lo cual utilizamos un bloque que llame al shell: %%Shell, con lo cual estaremos utilizando los comandos que se usan cuando se utiliza el shell

```shell
%%shell
#Seguimos el pipeline de FranciscoZorrilla Metagem, el cual es para establecer relaciones metabolicas en metagenomas,
#sin embargo, los primeros pasos para el filtrado de calidad se pueden seguir tal cual

#instalar y ejecutar miniconda

sudo apt-get update
wget https://repo.continuum.io/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
export PYTHONPATH="$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update --q conda
```

Para realizar la instalación de FastQC utilizamos el siguiente bloque 
```python
# Instalar FastQC para identificar problemas de calidad en los datasets. 
!conda install -c bioconda fastqc -y

```
Además podemos instalar fastp (puedes revisar la calidad ya sea utilizando fastqc o fastp, fastp también sirve para preprocesar en lugar de trimmomatic y cutapad)

```python
# Instalar fastp para identificar problemas de calidad en los datasets. 
!conda install -c bioconda fastp -y
```


Ahora instalamos también MultiQC, este agregará todos los análisis de calidad de los datos 
```python

# Instalar Multiqc que puede agregar y sumar todos los QC de los datos y alinear log data en un solo archivo 
!pip install multiqc
```
En caso de que querramos utilizar trimmomatic y no solo fastp podemos instalarlo a su vez.

```python
# Instalar Trimmomatic una herramienta flexible de trimming de lecturas para datos de Illumina NGS data 
! conda install -c bioconda trimmomatic -y
```

En caso de que queramos trabajar con transcriptomas (RNASeq) podemos utilizar Kallisto, usando las opciones:
Pseudobam

  --pseudobam outputs all pseudoalignments to a file pseudoalignments.bam in the output directory. This BAM file contains the pseudoalignments in BAM format,   ordered by reads so that each pseudoalignment of a read is adjacent in the BAM file.

  A detailed description of the SAM output is here.
 y GenomeBam

  --genomebam constructs the pseudoalignments to the transcriptome, but projects the transcript alignments to genome coordinates, resulting in split-read alignments. When the --genomebam option is supplied at GTF file must be given with the --gtf option. The GTF file, which can be plain text or gzipped, translates transcripts into genomic coordinates. We recommend downloading a the cdna FASTA files and GTF files from the same data source. The --chromosomes option can provide a length of the genomic chromosomes, this option is not neccessary, but gives a more consistent BAM header, some programs may require this for downstream analysis. kallisto does not require the genome sequence to do pseudoalignment, but downstream tools such as genome browsers will probably need it.
  
Es un programa para cuantificar la abundancia de &#x1F534; **Transcritos** a partir de datos de RNA-Seq en masa y unicelulares, o más generalmente de secuencias objetivo utilizando lecturas de secuenciación de alto rendimiento. Se basa en la novedosa idea de la pseudoalineación para determinar rápidamente la compatibilidad de las lecturas con los objetivos, sin necesidad de alineación.


Una vista más detallada de kallisto la podemos encontrar aquí: [Manual Kallisto](http://pachterlab.github.io/kallisto/manual.html)

```python
# Puede que Kallisto no sea necesario: es para cuantificar abundancias de transcritos de datos de RNA-seq, o 
#o de forma mas general secuencias blanco usando high-throughput sequencing reads. Podemos omitirlo en este cuaderno 
! conda install -c bioconda kallisto -y
```

Para utilizar metagem tenemos que clonar el repositorio del autor. Metagem nos permite reconstruir modelos metabólicos directamente de metagenomas. Se resumirá un poco más adelante.
[Artículo de Metagem](https://academic.oup.com/nar/article/49/21/e126/6382386)


```shell
%%shell
#clonar el repositorio de metagem y despues moverse al directorio

git clone https://github.com/franciscozorrilla/metaGEM.git
cd metaGEM

#establecer mamba y metagem
/root/miniconda/bin/conda create --quiet --prefix ./envs/mamba mamba -c conda-forge && source /root/miniconda/bin/activate envs/mamba
mamba env create --quiet --prefix ./envs/metagem -f envs/metaGEM_env.yml

#activar metagem en el ambiente conda e instalar las herramientas pip
source /root/miniconda/bin/activate envs/metagem
pip install --quiet --user memote carveme smetana 

#desactivar metagem y activar el ambiente mamba
source /root/miniconda/bin/deactivate && source /root/miniconda/bin/activate envs/mamba

#Establecer metawrap y prokka-roary
mamba env create --quiet --prefix ./envs/metawrap -f envs/metaWRAP_env.yml
mamba env create --quiet --prefix ./envs/prokkaroary -f envs/prokkaroary_env.yml 
```

## 2. Paquetería Python

En caso de que necesitemos utilizar paquetes para graficar y para trabajar con archivos csv debemos importar los paquetes matplotlib y pandas

```python
import matplotlib.pyplot as plt
import pandas as pd #llamar pandas
```

## 3. Montar Drive en colab
Para montar el drive utilizamos el siguiente bloque, esto es necesario para poder utilizar los archivos que se guarden en el drive y para guardar outputs en el mismo:

```python
#Montar el drive para trabajar con los archivos dentro del drive
from google.colab import drive
drive.mount('/content/drive')

```
Ya establecida toda la paquetería e instalados todos los pipelines podemos comenzar, teniendo en cuenta que nuestros archivos de secuencias (Pair-End o Single-End) deben estar guardados en nuestro drive.

## 4. Comencemos explorando las reads:

Las Lecturas comúnmente se encuentran en formato FASTQ, muy similar al FASTA solo que contiene caracteres alfanuméricos que dan indicio de calidad asociada con cada nucleótido. Una imagen de la estructura de un archivo fastq la podemos observar más adelante en el apartado del pre-procesamiento.

[Estructura FASTQ](#5-pre-procesamiento-analisis-de-calidad-usando-fastqc-y-fastp)

Son Archivos (de texto o ficheros) muy grandes que no se pueden leer. Ciertos comandos te permiten observar aspectos clave (Head/tail/more - trabajando con ficheros en Linux). 

Para conocer un poco nuestras secuencias podemos utilizar script de bash. En este caso trabajamos con los siguientes archivos:


![Nuestro drive y la carpeta conteniendo las lecturas:](https://user-images.githubusercontent.com/13104654/204870524-c62fcf00-b097-4758-b0a8-406bee184927.png)

Nótese que los archivos están comprimidos y se pueden seguir trabajando de esta manera.

### a) Explorando el contenido 

Utilizando el bash podemos observar encabezados y colas de archivos, con el comando zcat podemos trabajar con archivos comprimidos, en caso de no tenerlos comprimidos simplemente usamos cat. Como las lecturas son pair end, lo realizamos con ambos archivos.

```bash
#explorar el encabezado del archivo
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | head
```
```shell
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | head
```

```shell
#explorar el final para las lecturas R1 y R2
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | tail -n 4
```
```shell

%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | tail -n 4
```

### b) Revision del numero y longitud de secuencias 

> (También revisamos que se encuentren pareados, esta parte se puede omitir ya que fastqc te muestra algunos de estos resultados)

El siguiente bloque permite revisar si los archivos están pareados
```bash
##Para explorar que los archivos están pareados
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | wc -l
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | wc -l
```
si queremos saber cual es el número de secuencias usamos el siguiente (aunque es inespecífico), de nuevo usamos zgrep por que es un archivo comprimido.

```bash
##Explorar la cantidad de secuencias (debería ser 3949727, número de líneas/4)
%%bash
zgrep '^@' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | wc -l
zgrep '^@' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | wc -l
```
para hacer que el anterior bloque haga una búsqueda más específica tenemos

```bash
##Hay que ser más específicos para revisar la cantidad de secuencias
%%bash
zgrep '^@M02521' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | wc -l
zgrep '^@M02521' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | wc -l
```

Exploramos la longitud de secuencias. Para ello podemos usar awk, el cual es un comando (lenguaje de programación) que te permite trabajar operaciones matemáticas con ficheros de texto (filtrando). Más acerca de awk lo puedes encontrar [aquí](http://www.linuxfocus.org/English/September1999/article103.html#lfindex0) o [acá](https://geekland.eu/uso-del-comando-awk-en-linux-y-unix-con-ejemplos/).

Para cada línea de secuencia podemos contar cada caracter usando el parámetro NR (número de registros) y usando el contador.

```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}; print "total reads: " counter}'

```
y podemos añadir parámetros para imprimir en txt

```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}}' | sort -n | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/read_length1.txt
```

el txt generado lo podemos importar y transformar a csv para después graficarlo en matplotlib, para lo cual se pueden aplicar los siguientes bloques:

```python
read_file = pd.read_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_length1.txt',header=None)
read_file.to_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_lengthR1.csv', index=None,header=None)
```
Generar el df tipo decimal

```python
read_file[0] = read_file[0].astype(str)
df = read_file[0].apply(lambda x: x.strip(" "))
df = pd.DataFrame(df.str.split(' ',2).tolist(),
                                 columns = ['fips','length', 'ocurrences'])
df
df=df.astype(float)
```

y graficamos
```python
##Graficar la longitud de secuencia contra las ocurrencias
df.plot(x= 'ocurrences', y= 'length', kind ='line')
plt.show()
```
lo anterior lo aplicamos también en Reverse 

```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}; print "total reads: " counter}'
```

```bash
#Output de distribución de longitudes de secuencias y graficar en matplotib
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}}' | sort -n | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/read_length2.txt
```

```python
read_file2 = pd.read_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_length2.txt',header=None)
read_file2.to_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_lengthR2.csv', index=None,header=None)
```

```python
read_file2[0] = read_file2[0].astype(str)
df2 = read_file2[0].apply(lambda x: x.strip(" "))
df2 = pd.DataFrame(df2.str.split(' ',2).tolist(),
                                 columns = ['fips','length', 'ocurrences'])
df2
df2=df2.astype(float)
```

```python
##Graficar la longitud de secuencia contra las ocurrencias
graficaR2 = df2.plot(x= 'ocurrences', y= 'length', kind ='line')
plt.show()
plt.savefig("ReadlenghtR2.jpg", bbox_inches='tight')
```
Lo siguiente puede ser opcional.

### c) bloques para mostrar solo secuencias o secuencias repetidas

```bash
#Desplegar solo la información de la secuencia
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | head -n 10
```

```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | head -n 10
```
podemos considerar un output txt para las secuencias repetidas

```bash
#Observar barcodes-secuencias que se presentan más frecuentemente
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | sort | uniq -c | sort -n -r | head -n 10
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | sort | sort -n -r | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/Sec_rep1.txt
```
```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | sort | uniq -c | sort -n -r | head -n 10
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz | awk 'NR % 4 == 2 {print;}' | sort | sort -n -r | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/Sec_rep2.txt
```

## 5. Pre procesamiento: analisis de calidad usando FASTQC y FASTP

El archivo FastQ que habíamos abordado [anteriormente](#4-comencemos-explorando-las-reads), 
consta de cuatro líneas:

&#x1F535; 1. Nombre de la secuencia (header - id del secuenciador, coordenadas del spot, flow-cell, adaptador, etc.)

&#x1F535; 2. Secuencia

&#x1F535; 3. Espaciador (+), que opcionalmente puede contener más texto, y

&#x1F535; 4. Valores de calidad: Q Score - alfanumérico. Correspondientes a los nucleótidos en línea 2
    
Tal como se observa en la siguiente imagen:
![Estructura FastQ](https://user-images.githubusercontent.com/13104654/204933554-fef6cf9a-f8d4-4e52-ad1b-a831d5bfdd92.png)

Los scores de calidad (Q) representados en la cuarta línea, que corresponden a caracteres ASCII, reflejan, en escala logarítmica, la probabilidad de que esta base en particular fuera llamada incorrectamente (P<sub>error</sub>).

```math
Q = -10 log_{10} P               
```
```math
P_{error} = 10^{-Q/10}
```
Lo cual corresponde generalmente a:

![Phred Quality score](https://user-images.githubusercontent.com/13104654/205389945-0d25371a-f02b-4db0-8e80-07e2cedf0e41.png)

Hay que considerar además el método de codificado dependiendo de plataforma:

![Quality scores dependiendo de la plataforma de secuenciación](https://user-images.githubusercontent.com/13104654/205388609-2d6df438-0aea-4a9e-a612-37563d0e83e6.png)

Así, podemos además, aplicar otro código para determinarr si el score de nuestras lecturas corresponde a Phred+33, Phre+64 o Solexa+64

```bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz | head -n 10000 |\
  awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
  awk 'BEGIN{min=100;max=0;} \
      {for(i=1;i<=NF;i++) \
          {if($i>max) max=$i; \
               if($i<min) min=$i;}}END \
          {if(max<=74 && min<59) \
                     print "Phred+33"; \
           else \
           if(max>73 && min>=64) \
                     print "Phred+64"; \
           else \
           if(min>=59 && min<64 && max>73) \
                     print "Solexa+64"; else print "Unknown score encoding!";}' --- [Table of contents](/introduction/terminology_index.html)
```

Para entender un poco más acerca de [calidades](https://maq.sourceforge.net/qual.shtml) y los archivos fastq podemos ir a la documentación oficial de [FastQ](https://maq.sourceforge.net/fastq.shtml) y este [artículo](https://pubmed.ncbi.nlm.nih.gov/9521921/), además de esta [nota técnica](https://www.illumina.com/documents/products/technotes/technote_Q-Scores.pdf) de Illumina.


Tomando en cuenta lo anterior, podemos utilizar algunas herramientas para identificar problemas de calidad en las lecturas, con el fin no solo de mantener las secuencias adecuadas si no también de reducir el tamaño del archivo, evitar contaminación, etc.

Entre las herramientas a utilzar tenemos:

# 5.1 Fastqc
Para correr FastQC en los archivos de secuencias dentro de google colab usamos el siguiente bloque de código:

```python
# Pre-alignment QA 
!fastqc /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz
!fastqc /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz
```
Los archivos de salida son html y abarcan las siguientes evaluaciones:

## 5.1.1. Estadística simple

El primer apartado de estadística simple contiene el nombre del archivo, el número total de secuencias, si existen secuencias de mala calidad, longitud de las secuencias y el contenido de GC. Ya se había abordado en los scripts anteriores la determinación del total y longitud de secuencias, por lo que nos sirve para corroborar. En este apartado no se generan warnings o fails.


![Panorama general y vista de los estadísticos iniciales](https://user-images.githubusercontent.com/13104654/205708653-93a21bca-be14-44e7-839a-fba67d08786e.png)

## 5.1.2.  Calidad de secuencias

En este apartado se muestra una revisión del rango de los valores de calidad a través de todas las bases en cada posición en los archivos FASTQ. 
En cada posición se muestra una boxplot con bigotes. En la gráfica podemos definir una mediana (línea roja central), rangos intercuartiles 25-75%(cajas amarillas), los bigotes representan los puntos del 10 y 90% y la calidad media (línea azul). 

![calidad de secuencias por base para la lectura R1 de P69](https://user-images.githubusercontent.com/13104654/205745971-2a852136-b15d-431f-8720-de0edb5af83c.png)


El eje y corresponde a los *scores* de calidad, el cual es dividido con un fondo verde, naranja y rojo, siendo el fondo verde para los mejores *scores*, el naranja para los *scores* de no tan buena calidad y el rojo a los de calidad pobre (entre más alto el *score* mejor). 

Es normal para todas las plataformas que conforme avance la corrida la calidad disminuya. En esta parte se puede generar un warning, en el caso de que el cuartil de cualquier base sea menor de 10 o si la mediana es menor de 25 y fails si el cuartil es menor de 5 y la mediana menor de 20.

En el análisis que se realizó del genoma de P69 se muestra que la calidad media no cae al fondo rojo y su comportamiento es típico disminuyendo la calidad con la corrida. 
Se observa en la última parte solo la mediana de una caja en el umbral de 20, y los cuartiles menores de 10 pero mayores de 5, por lo tanto solo se lanza un warning, el cual será resuelto al realizar el preprocesamiento.  


## 5.1.3.  Calidad de secuencias por pozo (*flowcell*) 

En este apartado se muestra un heatmap de las pérdidas de calidad posicional, es decir se grafica la calidad de las llamadas de bases contra la posición física del secuenciador de la cual provienen. En los secuenciadores illumina, el área del *flowcell* se divide artificialmente en franjas (superficie superior e inferior) y 
estas se subdividen en mosaicos (áreas arbitrarias que se analizan por separado en la canalización de llamadas). Observar la calidad por mosaico identificará este tipo de errores. Se espera siempre la pérdida de calidad conforme los ciclos se incrementan, por ello resulta útil realizar una normalización para representar la calidad. Así un buen gráfico se observará en blanco (sólido azul rey brillante). 
De esta forma, los problemas se pueden manifestar de varias formas, en la prueba de P69:

![P69 calidad de secuencia por flowcell R1](https://user-images.githubusercontent.com/13104654/210288959-a6367307-2827-4ff1-b6ec-2dfb468564c0.png)

Puede haber una pérdida aleatoria de calidad en las posiciones y ciclos, lo cual indica un problema generalizado en el que más comúnmente se sobrecarga la celda. En el análisis de calidad de la imagen de la cepa P69 se observa algo de este problema generalizado con la corrida aunque algo menos intenso. Resulta un tanto problemático si mosaicos aparecen en áreas amplias y localizadas  del *flowcell*.

Si se pueden observar las pérdidas de calidad en mosaicos específicos entonces es posible removerlos del análisis *downstream*, lo que resultaría algo problemático de estar al inicio de la corrida. 

Ya que no sabemos cuantas lecturas son afectadas entonces la mitigación en este caso (P69) podría resultar un problema (al remover lecturas se podría perder información).

## 5.1.4.  Scores de calidad por secuencia 
Este apartado muestra un gráfico del número de secuencias (Y) contra la escala logarítmica del Phred (X), indicando la calidad de cada lectura:

![Calidad por secuencia](https://user-images.githubusercontent.com/13104654/210292722-b817acd1-1c04-415f-b470-91bc1f2aa7fd.png)

Phred score = 30 indica un rango de error de 1 base en 1000, o una exactitud de 99.9%.
Phred score = 40 indica un rango de error de 1 base en 10,000, o una exactitud de 99.99%.
Si el Phred score es < 27 se obtendrá un warning y por debajo de 20 se dará un fail. 

En el caso de P69, el promedio de calidad es 36, lo cual es bueno.

## 5.1.5.  Contenido de bases por secuencia

En el caso del contenido de bases por secuencia, este apartado nos muestra, como el nombre lo indica, la composición porcentual de las bases en cada posición de la secuencia. Como habría que esperar esta composición debe permanecer estable en todos los ciclos, claro tomando en cuenta que el contenido de bases puede variar dependiendo de ciertos factores como la especie. En algunos casos podemos observar un sesgo en las primeras partes de la corrida, como es el caso del presente análisis P69. Se observa claramente que dicho sesgo se disipa en el resto de la corrida. 


![Composición de secuencia](https://user-images.githubusercontent.com/13104654/210295270-29332fcd-2e65-4814-9439-5f7f743b6ab8.png)


En este apartado se puede generar una  :warning: **alerta**, si el contenido de bases varía más del 10% en cualquier posición, y generará un :x: **fail** si este porcentaje de variación es mayor al 20%.

La causa de este sesgo puede ser el paso de *priming* aleatorio en la producción de las librerías. A pesar de que los hexameros en el priming deben presentarse con igual frecuencia en el mix y deberían realizar el *prime* con eficiencia similar, en la realidad no se da el caso y ciertos hexámeros son favorecidos durante este paso.

¿Entonces, el sesgo tendría implicaciones en los análisis downstream?. Hay algunos puntos a tomar en cuenta:

- Es posible que como parte del sesgo haya un incremento en el *mis-priming* - ocasionando un número alto de *mis-called* bases al inicio de la secuencia, y, 
- Es posible que la selección del sesgo introducido tenga un efecto significativo en la capacidad de la librería de medir el contenido original debido a ciertas secuencias favorecidas.

Sin embargo estos puntos pueden no representar un gran problema ya que son fácilmente detectados,  algunos mencionan que pueden mitigarse por un *Trimming 5'*, sin embargo esto no es un arreglo. Ya que la composición sesgada es creada por la selección de fragmentos de secuenciado y no por errores de llamadas de bases, el único efecto del *trimming* es cambiar de tener una libreria que inicia en posiciones sesgadas a una que inicia más allá de dichas posiciones. 

La única forma de resolver este problema sería introducir nuevos kits de preparación de librerías con una menor disposición al sesgo en el paso del *priming*, sin embargo, a pesar de la advertencia,  no parece que haya consecuencias serias para los análisis posteriores, irónicamente en RNA-seq son más sospechosas las librerías que no presentan este artefacto.

En este [artículo](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085583) se documenta el *mis-priming* en RNA-seq

Algunas de las razones más comunes de un warning o fail en este apartado son:
Secuencias sobrerepresentadas, sesgo en la fragmentación y composición sesgada de las librerías (que a veces ocurre naturalmente).

En el presente caso también se observa una desviación al final, si se está analizando una biblioteca que ha sido recortada agresivamente por el adaptador, naturalmente introducirá un sesgo de composición al final de las lecturas, ya que las secuencias que coinciden con tramos cortos del adaptador se eliminan, dejando solo las secuencias que no coinciden. Por lo tanto, es probable que las desviaciones repentinas en la composición al final de las bibliotecas que han sufrido un recorte agresivo sean falsas.


## 5.1.6.  Contenido de GC por secuencia

Este apartado muestra en un plot, el contenido porcentual total de GC para todas las lecturas (número total de reads vs porcentaje de bases G y C por lectura), comparando contra una "distribución teórica" de GC's, asumiendo un contenido uniforme para todas las lecturas, el pico central corresponde al contenido de GC promedio del genoma subyacente. Dado que el contenido de GC del genoma no se conoce, el contenido modal de GC es calculado de datos observados y usado para construir la distribución de referencia.

![Contenido de GC](https://user-images.githubusercontent.com/13104654/210295287-8ed4d1b0-051f-498d-b5f9-d2cdbb5d60c3.png)

 Se observa un ⚠️**warning** si el 15% total de las secuencias caen fuera de la distribución normal.
 
 Se obtendrá un :x:**fail** si más del 20% (el manual FastQC indica 30%) de las secuencias están fuera de la distribución normal.
 Los fails son generalmente debidos a contaminación, frecuentemente por secuencias de adaptadores.

Una distribución de forma inusual podría indicar una librería contaminada o alguna otra clase de subset sesgado. Una distribución normal cambiada índica algún sesgo sistemático, el cual es independiente de la posición de la base. Si existe un error sistemático, este no será marcado como error por que no se sabe cual debería ser el contenido de GC del genoma.

Existen otras situaciones en las cuales una distribución inusual se puede presentar. Por ejemplo, con RNA seq puede haber una distribución mayor o menor del contenido medio de GC entre los transcritos, causando que el gráfico observado sea más amplio o más estrecho que una distribución normal ideal.

## 5.1.7.  Contenido de N por base 

Si un secuenciador no puede llamar una base con confianza suficiente entonces será sustituido normalmente con un N más que una llamada de base convencional.

Este módulo gráfica el porcentaje de "llamadas" de base en cada posición para las cuales una N fue considerada.

![Contenido de N P69](https://user-images.githubusercontent.com/13104654/210295307-9831c8c6-f696-4640-ad40-0efc91a27528.png)

Idealmente el contenido de N por base sería una línea plana en 0% sobre el eje Y, indicando que todas las bases han sido "llamadas".

  - Se recibe un :warning: **warning** si el contenido de N es igual o mayor de 5%,
  - Tendremos un :x: **fail** si el contenido de N es igual o mayor a 20%.

El análisis de R1 para P69 muestra el resultado ideal para este módulo.

## 5.1.8.  Distribución de la longitud de secuencia

Este gráfico muestra, tal como hicimos en [apartados anteriores](#b-revision-del-numero-y-longitud-de-secuencias), la distribución de los tamaños de fragmentos en el archivo analizado. En muchos casos esto solo produce un gráfico simple mostrando solo un pico de un solo tamaño, pero para archivos FASTQ con longitud variable, mostrara las cantidades relativas de cada tamaño de fragmento de secuencia. En este caso nuestro archivo R1 P69 muestra longitudes variables más pequeñas (30-269pb) que el pico de 299pb. 

![Longitud de secuencias R1 P69](https://user-images.githubusercontent.com/13104654/210295322-a0759d95-e0a6-42cc-b699-343741a5974e.png)


Algunos secuenciadores (y kits de secuenciación) generan fragmentos de longitudes ampliamente variables, otros pueden generar fragmentos de longitud uniforme.
Incluso en librerías con longitud uniforme, algunos *pipelines* cortarán secuencias para remover llamadas de bases de baja calidad del final o las primeras n bases si coninciden las primeras n bases del adapatador arriba del 90% (por defecto), algunas veces con n=1. 
Para secuenciación con Illumina, cada lectura debería ser del mismo tamaño (?). 

Este módulo arrojará un :warning:**warning** si hay cualquier variación en la longitud de las secuencias, el cual puede ser ignorado si se sabe que es normal para el tipo de datos que se tiene. 

Un :x: **fail** en este módulo significa que al menos una secuencia tiene longitud de 0. 
El análisis de R1 P69 obtiene un warning ya que hay una gran variabilidad en la longitud de las secuencias, lo cual puede cambiar al realizar el trimming.


## 5.1.9.  Niveles de duplicación de secuencias

En este módulo se grafican los niveles de duplicación de secuencias (eje x) contra el porcentaje de secuencias que muestran ese nivel de duplicación (eje y), y 

Hay dos líneas en el gráfico:
🔴 línea roja: Distribución para las secuencia de-duplicadas con las proporciones del conjunto de-duplicado las cuales provienen de diferentes niveles de duplicación en los datos originales.

🔵 línea azul: Distribución de los niveles de duplicación para en conjunto completo de secuencias. 

![niveles de duplicación de secuencias en R1 P69](https://user-images.githubusercontent.com/13104654/210295333-f2d9ca7d-0193-4483-b601-e077f178b6c3.png)

La gráfica de los niveles de duplicación de secuencias muestran en el eje x, el número de veces que una secuencia está duplicada, y en el eje y el porcentaje de secuencias que muestran ese nivel de duplicación. Normalmente un genoma tendrá un nivel de duplicación de secuencias de 1 a 3 para la mayoría de las secuencias, con sólo un puñado de lecturas teniendo un nivel más alto que este; la línea debería tener la forma inversa a una gráfica log.

En el presente análisis de R1 P69 se no se observan picos a la derecha de la gráfica y solo un bajo nivel de duplicación al inicio

Un alto porcentaje de duplicación de secuencias es un indicativo de contaminación.

Este módulo nos arrojará un :warning: **warning** si más del 20% de las secuencias son duplicadas.

Tendrémos un :x: **fail** si más del 50% de las secuencias están duplicadas. 

Un warning o fail pueden ser resultado de artefactos de PCR.

#### Más acerca de la duplicación:

>En una librería diversa la mayoría de las secuencias se presentarán solo una vez en el set final, un bajo nivel de duplicación puede indicar un muy alto nivel de coverage de la secuencia blanco, pero un alto nivel puede indicar una clase de sesgo por enriquecimiento ( por ejemplo en la amplificación por PCR).
Este módulo cuenta el grado de duplicación para cada secuencia en el conjunto y crea un plot mostrando el numero relativo de secuencias con diferentes grados de duplicación.

>Con el fin de reducir los requerimientos de memoria para este módulo, solamente las secuencias que se presentan en las primeras 200 000 en cada archivo son analizadas, pero esto debería bastar para obtener una impresión para los niveles de duplicación del archivo completo. 
Cada secuencia es rastreada al final del archivo para dar un conteo representativo del promedio del nivel de duplocación. 
Para reducir la cantidad de información en el gráfico final, cualquier secuencia con >10 duplicados son colocadas en esta categoría, por lo que no es inusual observar un leve incremento en esta categoría final. Si hay un gran incremento, significa que se tiene un alto número de secuencias con alto nivel de duplicación. 

>Debido a que la detección de la duplicación requiere de una coincidencia exacta de secuencias sobre la longitud completa de la secuencia, cualquier lectura por encima de 75pb de longitud son truncadas a 50pb para propósitos del análisis, aún así, lecturas más largas son más propensas a contener errores de secuenciamiento por lo cual incrementará artificialmente la diversidad observada y tenderá a subrepresentar las secuencias altamente duplicadas.

>Para datos del *Whole Genome Shotgun* se espera que cerca del 100% de las lecturas sean únicas (una sola vez en los datos de secuencia). La mayoría de las secuencias deberían caer hacia la izquierda del gráfico para ambas líneas. Esto indica una librería altamente diversa que no esta sobre secuenciada. Si la profundidad del secuenciamento es extremadamente alta (p. ej. >100x el tamaño del genoma) es inevitable que aparezcan duplicaciones de secuencias: en teoría solo hay un número finito de lecturas de secuencia completamente únicas las cuales pueden ser obtenidas de cualquier muestra de DNA ingresada.

>Subconjuntos de enriquecimiento más específicos, o la presencia de contaminantes de baja complejidad tenderán a producir picos hacia la derecha del gráfico. Estos picos de altos niveles de duplicación aparecerán más frecuentemente en la línea azul ya que conforman una mayor proporción de la librería original, pero usualmente desaparecen en el trazo rojo, ya que consiste de una porporción no significante del conjunto deduplicado. Si los picos persisten en la línea roja, entonces esto sugiere que hay un alto número de secuencias diferentes altamente duplicado lo que podría indicar ya sea un conjunto de contaminantes o una duplicación técnica severa.

>Es usualmente el caso para RNA seq donde existen algunos transcritos altamente abundantes y algunos con baja abundancia. Se espera que las lecturas duplicadas sean observadas para los transcritos de alta abundancia.


## 5.1.10. Secuencias sobre representadas

En el caso de este módulo:
- Si se calcula que alguna secuencia representa más del 0.1 % del genoma completo será etiquetada como una secuencia sobre-representada y se obtendrá un :warning: **warning**
- La presencia de secuencias que representan más del 1% del genoma dará como resultado un :x: **fail**.

![Sobrerrepresentación R1 P69](https://user-images.githubusercontent.com/13104654/210641941-b7fb8d5a-2bce-4183-afa8-bf31b0cf0096.png)

En el presente análisis no se presentaron secuencias sobre-representadas.


> Una librería normal contendrá un conjunto diverso de secuencias, ninguna de las cuales individualmente hace una fracción del completo. Encontrar que una sola secuencia se encuentra sobre representada en el conjunto o significa que es altamente significativa biológicamente, que la librería está contaminada o bien que no es tan diversa como se esperaba.

> FastQC enlista todas las secuencias que hacen más del 0.1% del total y por cada secuencia busca coincidencias en una base de datos de contaminantes comunes y reportará el mejor *Hit*. Los *Hits* deben ser de al menos 20pb en longitud y tener máximo un *mismatch*. Encontar uno no necesariamente significa que sea la fuente de contaminación pero puede apuntar en la dirección correcta. Muchas secuencias de adapadores son muy similares entre sí, por lo que podría tenerse una coincidencia técnicamente incorrecta.

> Los datos de RNAseq pueden tener algunos transcritos que son tan abundantes que se registran como secuencias sobre-representadas. 
Con los datos de DNA seq, ninguna secuencia debería presentarse con suficientemente alta frecuencia para ser listada, pero algunas ocasiones podemos encontrar un pequeño porcentaje de lecturas de adaptadores.

> Podemos hacer BLAST de la secuencia sobre representada, si [Blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) no nos proporciona respuesta, podemos utilizar [VecScreen](https://www.ncbi.nlm.nih.gov/tools/vecscreen/).


## 5.1.11. Contenido de adaptadores 

Este módulo busca secuencias específicas de adaptadores.

   - Una secuencia que representa más del 5% del total causará un :warning: **warning** en este módulo.
   - Una secuencia que represente más del 10% del total causará un :x: **fail**.

![Contenido de adaptadores R1 P69](https://user-images.githubusercontent.com/13104654/210295352-5c134059-dc4a-4e72-bee1-18e3ee3eadbb.png)

Nuestro análisis no muestra contaminación con secuencias de adaptadores, lo cual es ideal. 
Si existiera un número significativo de secuencias de adaptadores, se debe utilizar un programa para recortarlos y realizar el análisis de calidad nuevamente.

Otros gráficos relacionados pueden consultarse en [Documentación Contenido de K-mer FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html) y más problemáticas en [![QC Fail]](https://sequencing.qcfail.com/).


Lo anterior hay que correrlo para las lecturas R2. Luego podemos utilizar Multiqc para tomar ambos resultados de calidad o bien podemos usar fastp para realizarlo en un solo paso.

## 5.1.12 Multiqc

El siguiente bloque de código nos sirve para correr multiqc en colab.

```python
# Pre-alignment Multiqc summary file 
!multiqc .

import multiqc
#El análisis se hace sobre los archivos de fastqc para que te entregue un reporte conteniendo todo
multiqc.run('/content/drive/MyDrive/códigos/Secuencias_UADY/50-3_S19_L001_R2_001_fastqc.zip')

```


# 5.2 Fastp
Fastp es una herramienta que realiza el preprocesamiento y filtrado de calidad de forma paralela y soporta lecturas Single end y Paired end.
Más información se puede encontrar en el [repositorio](https://github.com/OpenGene/fastp#simple-usage) de los desarrolladores.

A comparación de FASTQC, fastp ofrece resultados tanto para los datos de prefiltrado como para los datos de post-filtrado, permitiendo una evaluación del efecto del filtro comparando directamente las gráficas y reporta sus resultados tanto en formato HTML como en formato JSON, siendo este último manualmente optimizado para facilitar su lectura (más acerca de la descripción en el [artículo](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)).
 
Para correr Fastp en los archivos de secuencias (con los datos para filtrado por defecto) dentro de google colab usamos el siguiente bloque de código:

```python
# Control de calidad y reporte 
!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz

```
Si se requiere establecer un límite de longitud para filtrado se utiliza -l, para establecer el nombre de los archivos de salida -j -h, más opciones [aquí](https://github.com/OpenGene/fastp#all-options)

```python
# Control de calidad y reporte 
!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz -R content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/fastp_report

#o también de esta forma
# Control de calidad y reporte 
!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz --json="HA1AB3SS04_S4_L1.json" --html="HA1AB3SS04_S4_L1.html" -l 150 --detect_adapter_for_pe -c --cut_right --cut_front -p --failed_out="failed_seqsPR69.fastq.gz"
```
En este caso Fastp puede detectar adaptadores y cortarlos lo que nos ahora tiempo, sin embargo esto se puede realizar con un script aparte utilizando trimmomatic 
para el filtrado de calidad, en este caso habría que correr nuevamente los análisis de calidad con Fastqc para ver como quedaron las secuencias.


#5.3 Trimmomatic

Para el filtrado y corte de adaptadores con trimommatic, se puede usar el siguiente bloque de código 
```python
# Trimming 
!trimmomatic PE -phred33 R1.fastq R2.fastq R1\_paired.fq.gz R1\_unpaired.fq.gz R2\_paired.fq.gz R2\_unpaired.fq.gz ILLUMINACLIP:contams\_forward\_rev.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

o bien podemos señalarle la secuencia de adaptadores. En el caso de PR69, las secuencias no tienen adaptador (o eso creemos) y tienen la misma longitud. 
En el presente caso no conozco la secuencia de adaptadores de la plataforma pero si se requiere, aqui podemos encontrar algunos adaptadores [https://github.com/ATGenomics/adapters](https://github.com/ATGenomics/adapters)

El siguiente bloque de código se aplica en una máquina local y no en google colab como hemos utilizado hasta ahora, podría funcionar si llamamos el shell y le indicamos alguna ruta específica para clonar las secuencias del repositorio. 

*Para secuencias Paired End:

```shell
cd $HOME/ git clone https://github.com/ATGenomics/adapters.git $HOME/bin/adapters

adapters="$HOME/bin/adapters/NexteraPE-PE.fa"
```

adaptándolo para colab podría ser de la siguiente manera:

```shell
%%shell

###Para no poner toda la dirección asignamos una variable, probar si esto funciona
#### Asignar nombres de archivos de lecturas
#fwd="/content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz"
#rev="/content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz"
#name="PR69"

#Primero creamos el directorio adaptadores
mkdir -p /content/drive/MyDrive/Analisis_Posdoc/PR69/adapters

#luego clonamos el repositorio dentro de esta carpeta
cd /content/drive/MyDrive/Analisis_Posdoc/PR69/adapters git clone https://github.com/ATGenomics/adapters.git $HOME/bin/adapters

#una vez clonado se define como una variable a la secuencia de adaptadores que corresponda la plataforma que estemos usando
adapters="/content/drive/MyDrive/Analisis_Posdoc/PR69/adapters/NexteraPE-PE.fa"
```

En una máquina local se activa el entorno qc 
`source activate qc`

Sin embargo en google colab esto no parece ser necesario puesto que simplemente se llama el shell o bash con la correspondiente función


```shell
%%shell
trimmomatic PE -phred33 -threads 16 \ fwd \ rev  \ /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.trim.fq.gz \ /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.trim.fastq.gz /Analisis_Posdoc/PR69/HA1AB3SS04_S4.trim.fq.gz \ ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:90 CROP:150
```
Es muy importante señalar que en este caso hay que tomar en cuenta el orden en el que se colocan los parámetros.

Trimmomatic realiza un *trimming* de calidad adaptativo, cortado de cabeza y cola y remoción de adaptadores. Se puede revisar la documentación y bajar el programa [aquí](http://www.usadellab.org/cms/index.php?page=trimmomatic).

Una de las ventajas del programa es que permite trabajar con secuencias *Paired end*, reteniendo solamente pares coincidentes.
Otra ventaja es que permite coincidencias parciales y *overlapping* para la búsqueda de adaptadores.

Las opciones que podemos utilizar son las siguientes:



## 5.3.1 Eficiencia y formato

>Las siguientes se usan siempre antes de la invocación de los archivos de entrada y salida

- *threads*: este ajuste modifica el número de "hilos" de CPU que Trimmomatic debería usar en las computaciones. Una computadora típicamente tiene cerca de 2 núcleos. los cuales deberían corresponder a una cantidad de 4 hilos disponibles. 
- *phred*:  [-phred33 	-phred64]:  Este ajuste le dice al programa que codifique el archivo

> A partir de aquí son las opciones que van después de la invocación de *Inputs/Outputs* (estas opciones se presentan en mayúsculas) 

Opciones para cambiar la codificación (ver en [Apartado 5](#5-pre-procesamiento-analisis-de-calidad-usando-fastqc-y-fastp)):
Si se requiere leer la codificación de un tipo y sacar la codificación de uno diferente, éstas opciones son las que se necesitan utilizar.

- TOPHRED33: Convierte *scores* de calidad a Phred-33
- TOPHRED64: Convierte *scores* de calidad a Phred-64 

### 5.3.1.1 Cortado (*Cropping*)

Trimmomatic cuenta con varias opciones que pueden ser usadas simultáneamente o no:

-LEADING: Corta bases del inicio de una lectura, si está por debajo del umbral de calidad - adaptativa
-TRAILING: Corta bases del final de una lectura, si está por debajo del umbral de calidad - adaptativa
-CROP: Corta la lectura a una longitud específica
-HEADCROP: Corta el número específicado de bases del inicio de una lectura.

LEADING y TRAILING son cortado adaptativo, lo que significa que cortarán el inicio/fin de las lecturas si fallan la calidad especificada. Lo anterior difiere de CROP y HEADCROP, los cuales podrían cortar a una longitud o número de bases específicas (respectivamente), en este caso, el programa realizara el corte para todas las lecturas.

-MINLEN: se deshará de todas las lecturas que caen bajo una longitud especificada.

### 5.3.1.2 *Trimming* de calidad adaptativo

-SLIDINGWINDOW: realiza un trimming en una ventana de deslizamiento, cortando una vez la calidad promedio dentro de u
    SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.

It takes two values like SLIDINGWINDOW:4:15 which means “Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15”
Adapter trimming

Finally, trimmomatic will take a file with the sequences of your adapters and will trimm them out. It follows the following call: ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>. From their docs:

        fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used. See the section below or use one of the provided adapter files
        seedMismatches: specifies the maximum mismatch count which will still allow a full match to be performed.
        palindromeClipThreshold: specifies how accurate the match between the two ‘adapter ligated’ reads must be for PE palindrome read alignment.
        simpleClipThreshold: specifies how accurate the match between any adapter etc. sequence must be against a read.



ILLUMINACLIP:${adapters} eliminamos adaptadores con cierta frecuencia

![trimmomatic_adapter](https://user-images.githubusercontent.com/13104654/211375856-b34becba-e0e4-450d-8d0b-06552b13b296.png)


 SLIDINGWINDOW: 4:20 cuatro nucleótidos en promedio tienen una calidad menor a 20 se elimina la secuencia (incluyendo el par).
 MINLEN: mínimo de largo 130.
 CROP: Cortar la lectura a una longitud determinada



## 6. Ensamble *De Novo*

















