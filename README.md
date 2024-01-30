# MINERÍA GENÓMICA 
## CALIDAD, ENSAMBLE Y ANOTACIÓN

Contiene pipeline para ensamble y minería genómica de shotgun sequence con google colab

El protocolo de manera general:

![Diagrama de flujo_Análisis](https://user-images.githubusercontent.com/13104654/211899391-66c4a856-e193-44f2-baae-0ad84895b78a.png)

![PipelineA](https://user-images.githubusercontent.com/13104654/211914431-9c9197e0-58b6-4e54-a068-221093935a43.png)

![PipelineB](https://user-images.githubusercontent.com/13104654/211915131-852fbab8-fc34-4ba3-a3f7-d91135098e47.png)

Antes que nada tenemos que instalar las herramientas que usaremos en la nube.
Desde el cuaderno establecido primeramente se instalan todos los paquetes que se usarán y al final se monta el drive en el que se estará trabajando. Es preferible que esto se haga desde el inicio ya que cuando se instala un nuevo paquete se reinicia el entorno y lo que anteriormente llamamos ya no estará disponible (cada vez que se quiera hacer el procedimiento ya que el cuaderno/entorno se ha cerrado hay que instalar todo de nuevo y hacer el montaje del drive). 

## 1. Instalacion de Herramientas

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
## Instalacion de BWA y SAMtools

```python

!conda install -c bioconda bwa -y
!conda install -c bioconda samtools -y

```

## Instalacion de Megahit y SPAdes
```python
#Para instalar megahit
!conda install -c bioconda megahit

#para instalar SPAdes
!conda install -c bioconda spades=3.15.5 -y
```

## Instalacion de Quast

```python
#Para instalar Quast
!conda install -c bioconda quast
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
##Explorar la cantidad de secuencias (debería ser #####, número de líneas/4)
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
%%bash
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

## 5.1.2. Calidad de secuencias por base

Este es el valor de confianza de base con base en el Phred score que designa las series de score de calidad de las bases completas en su respectica locación en el archivo. un valor más allá de Q30 es considerado bueno, mientras que uno arriba de Q20 es generalmente aceptado.
En este apartado se muestra una revisión del rango de los valores de calidad a través de todas las bases en cada posición en los archivos FASTQ. 
En cada posición se muestra una boxplot con bigotes. En la gráfica podemos definir una mediana (línea roja central), rangos intercuartiles 25-75%(cajas amarillas), los bigotes representan los puntos del 10 y 90% y la calidad media (línea azul). 

![calidad de secuencias por base para la lectura R1 de P69](https://user-images.githubusercontent.com/13104654/205745971-2a852136-b15d-431f-8720-de0edb5af83c.png)


El eje y corresponde a los *scores* de calidad, el cual es dividido con un fondo verde, naranja y rojo, siendo el fondo verde para los mejores *scores*, el naranja para los *scores* de no tan buena calidad y el rojo a los de calidad pobre (entre más alto el *score* mejor). 

Es normal para todas las plataformas que conforme avance la corrida la calidad disminuya. En esta parte se puede generar un warning, en el caso de que el cuartil de cualquier base sea menor de 10 o si la mediana es menor de 25 y fails si el cuartil es menor de 5 y la mediana menor de 20.

En el análisis que se realizó del genoma de P69 se muestra que la calidad media no cae al fondo rojo y su comportamiento es típico disminuyendo la calidad con la corrida. 
Se observa en la última parte solo la mediana de una caja en el umbral de 20, y los cuartiles menores de 10 pero mayores de 5, por lo tanto solo se lanza un warning, el cual será resuelto al realizar el preprocesamiento.  


## 5.1.3.  Calidad de secuencias por pozo (*flowcell*) 

![illumina_flowcell](https://user-images.githubusercontent.com/13104654/212420432-e0c77336-d557-4a28-a25e-6e03ef1ab5a8.png)

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
multiqc.run('/content/drive/MyDrive/códigos/*_fastqc.zip') #Recordar actualizar la dirección de donde se encuentren los archivos de salida de qc

```
Para observar el html es con el siguiente bloque:
```python
import IPython

# Best if using Google Colab
IPython.display.HTML(filename='./multiqc_report.html')
#IPython.display.HTML(filename='/content/drive/MyDrive/códigos/Secuencias_UADY/multiqc_report.html') #Para dirección específica, marcando error de momento 

```

Es una herramienta que busca todos los archivos de control de calidad de nuestras secuencias y las resume en un solo reporte, además nos permite subrayar, por ejemplo,
ciertas muestras dentro del reporte, cambiar nombres o bien bajar las gráficas (interactivas) generadas a diferentes resoluciones. 
(es mejor, de momento, tratar de correr esta herramienta en entorno local (de preferencia en linux, unix o wsl) por que no he averiguado como correrla desde colab y que no la coloque en la ventana del cuaderno de trabajo si no que entregue el output externo  ❗📥 ☺️ ***Actualización: el código funciona adecuadamente y el html también se guarda en el entorno del drive.***).
 
# 5.2 Fastp
Fastp es una herramienta que realiza el preprocesamiento y filtrado de calidad de forma paralela y soporta lecturas *Single end* y *Paired end*.
Más información se puede encontrar en el [repositorio](https://github.com/OpenGene/fastp#simple-usage) de los desarrolladores.

A comparación de FASTQC, fastp ofrece resultados tanto para los datos de prefiltrado como para los datos de post-filtrado, permitiendo una evaluación del efecto del filtro comparando directamente las gráficas y reporta sus resultados tanto en formato HTML como en formato JSON, siendo este último manualmente optimizado para facilitar su lectura (más acerca de la descripción en el [artículo](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)).
 
Para correr Fastp en los archivos de secuencias (con los datos para filtrado por defecto) dentro de google colab usamos el siguiente bloque de código:

```python
# Control de calidad y reporte 
!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz

```

## 5.2.1 Summary

En la salida primeramente podemos observar un resumen general antes y después del filtrado, así como el resultado del filtro
aplicado, esto lo muestra para los dos archivos de lecturas:

![fastp_summary](https://user-images.githubusercontent.com/13104654/213575454-29eaae55-82c1-4637-b1a5-895cd369057a.png)

En el caso de PR69, con los filtros por defecto, podemos observar que las lecturas de baja calidad contaron para un 3.6% del total, 
sin un porcentaje significativo de N ni lecturas cortas.

## 5.2.2 Adaptadores o mal ligado 

La siguiente sección muestra las ocurrencias de adaptadores de **ambos** archivos de lecturas.

![adapters_fastp](https://user-images.githubusercontent.com/13104654/213586493-70e52e97-87c2-465c-b510-8b91969bd51b.png)

Para PR69 muestra un bajo porcentaje de adaptadores (0.04% R1 y 0.4% R2). 

En este caso Fastp puede detectar adaptadores y cortarlos lo que nos ahorra tiempo, sin embargo esto se puede realizar con un script aparte utilizando [Trimmomatic](#62-trimmomatic) (o con [cutadapt](#631-cutadapt)) 
para el filtrado de calidad, en este caso, habría que correr nuevamente los análisis de calidad con Fastqc para ver como quedaron las secuencias.

## 5.2.3 Estimación del tamaño de inserto

En este apartado se muestra la distribución del porcentaje de lecturas (eje y) contra el tamaño de las lecturas (eje x) en un gráfico interactivo,
podemos modificar el tamaño de los ejes y hacer zoom.
Esta estimación toma en cuenta el overlap de las lecturas *Paired end*.

En el caso de PR69 se encuentran 56.33% de lecturas no sobrelapadas por lo que podrían ser de tamaño <30 o >572 o bien con gran cantidad de errores de secuenciación.

![Insert Size Distribution](https://user-images.githubusercontent.com/13104654/213606759-021f3824-8b43-470b-917c-17cecc3d64d9.png)

## 5.2.4 Antes del Filtrado
Las secciones siguientes muestran gráficos interactivos, antes del filtrado de calidad,  correspondiendo primero a la **calidad**, 
la gráfica es equivalente a la mostrada por fastqc en la sección 5.1.2 [Calidad de secuencias por base](#512-calidad-de-secuencias-por-base), 

![Calidad_R1](https://user-images.githubusercontent.com/13104654/213797876-82e8a084-a316-42d6-9d73-8d938506a170.png)


luego se muestra el gráfico de los índices del contenido de bases contra la posición, incluido el contenido de Ns

![Contenido de bases_R1](https://user-images.githubusercontent.com/13104654/213798468-0381814c-c3a1-4790-ab4d-3a4e5edf345e.png)

 y finalmente un heatmap con el conteo de K-meros, 
donde las áreas más oscuras representan cuentas mayores. Lo anterior primero para R1 y después para R2.

![Conteo de Kmer](https://user-images.githubusercontent.com/13104654/213798539-eacad1b7-1c3e-4d40-88aa-85b0d4932ea2.png)

La idea de los *K meros* es simple, se crea una ventana de longitud *k* y se desliza tomando un caracter al tiempo. Si la longitud de una secuencia de DNA dada es N, entonces tendremos:

```math
Total K-mer = N - k + 1               
```
Usualmente se buscan tres tipos de frecuencias, la cuenta total (que tantas veces aparece un *k-mero* en una secuencia dada), la cuenta diferida (si ha aparecido o no sin importar cuantas veces) y la cuenta única (aquellas que solo han aparecido una vez). 
El conteo de *k-meros* es útil para el ensamble, clustering y alineamientos, así como para la corrección de errores en secuenciado, la estimación del tamaño del genoma y la identificación de repeticiones. Computacionalmente es un proceso exigente.

![k-mer](https://user-images.githubusercontent.com/13104654/213822121-8dacf7ca-1e63-4d6d-9094-e6a807141c37.png)


## 5.2.5 Después del Filtrado
Luego aparece una sección con gráficos para después del filtrado (que puede ser por defecto o definido de acuerdo a lo que se observe en los datos).

El html se encuentra en la carpeta de Drive que se indicó en *colab*, y la salida de archivos ya filtrados también se encontraran donde se indicó, con este archivo fastq podemos continuar con los siguientes análisis *Downstream*. 

fastp también cuenta con una *flag* para realizar el *merge* de las dos lecturas. Esto una vez que se cuenta con las secuencias ya filtradas.

# 6. Pre-procesamiento: Filtrado de calidad 

Ahora que se tiene conocimiento acerca de los datos crudos, es importante usar estaa información para limpiar y *Trimmear* las lecturas para mejorar la calidad general antes del ensamble. Hay cierto número de herramientas disponibles para esta tarea (a varios grados), pero necesitamos lidiar con lecturas pareadas (en caso de tener lecturas *paired end*, como es el presente caso). Si uno de los *ends* de un par es removido, la lectura huérfana necesita colocarse en un archivo separado de "Lecturas huérfanas", lo cual mantiene el orden de pareado de las lecturas en los archivos para que el programa de ensamble las pueda usar correctamente.
Entre las herramientas más comunes disponibles son Trimmomatic, cutadapt, PRINSEQ, QC Chain entre otras, esta tarea también la puede realizar fastp modificando las flags correspondientes.

# 6.1 Fastp (filtrado)

En el caso del filtrado utilizando fastp podemos realizar lo siguiente.
Si se requiere establecer un límite de longitud para filtrado se utiliza -l, para establecer el nombre de los archivos de salida -j -h, más opciones [aquí](https://github.com/OpenGene/fastp#all-options)

### Longitud mínima de lectura
El valor más apropiado para este parámetro dependerá de los resultados del reporte de FastQC/Fastp, específicamente la longitud de alta calidad en el gráfico de la sección [Per Base Sequence Quality](#512-calidad-de-secuencias-por-base) y segunda y quintas gráficas en la sección de [Antes del filtrado](#524-antes-del-filtrado) de Fastp.

De nuestros resultados, podemos establecer este mínimo en 36

```python
# Control de calidad y reporte 
#!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz -R content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/fastp_report -a --detect_adapter_for_pe

#o también de esta forma
# Preprocesamiento 
!fastp -i /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz  -O /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz --json="HA1AB3SS04_S4_L1.json" --html="HA1AB3SS04_S4_L1.html" -l 36 --cut_right --cut_front -c -m --merged_out /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_mer.fastq.gz --unpaired1 /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_up1.fastq.gz --unpaired2 /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_up2.fastq.gz --failed_out /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_fout.fastq.gz
```
> Podemos usar --detect_adapter_for_pe antes de -- cut right por si hay adaptadores 

> Cut_right equivale a SLIDING WINDOW en Trimmomatic mueve una ventana deslizante desde el frente al final, si encuentra una ventana con calidad media < *Threshold*, se deshace de las bases en la ventana y la parte derecha, luego para, este valor lo podemos establecer a 4:20 (4pb/20 Phred) en Trimmomatic, aquí dichos valores se establecen por defecto, a menos que especifiquemos lo contrario con otros flags (--cut_right_window_size y --cut_right_mean_quality). 

> Cut_front se mueve del frente 5' a la cola, cortando las bases en la ventana que no alcanzan la calidad media.

> -c activa la corrección de bases en regiones sobrelapadas (solo para lecturas PE). 

> -m para inputs *paired end*, combina cada par de lecturas en una sola si se encuentran sobrelapadas. Las lecturas que se combinan serán escritas en el archivo dado por --merged_out, las lecturas sin combinar se especifican en --out1 y --out2. El modo combinado por defecto se encuentra desactivado.

> En el caso de unpaired1 y unpaired2 para PE, si la lectura1 pasa QC pero la lectura2 no será escrita en unpaired1, viceversa para unpaired2. Si unpaired1 y unpaired2 son la misma, ambas serán escritas en el mismo archivo.

>  en el flag de Failed_out se escriben las lecturas que no pasan los filtros.

CH606 Fastp bloque

```python

!fastp -i fwd.fq.gz -I Reverse.fq.gz -o fwdR1_filtrada.fq.gz  -O RvR2_filtrada.fastq.gz -l 100 --cut_right --cut_front -c -D --dedup -m --merged_out salida_merge.fq.gz --unpaired1 up1.fq.gz --unpaired2 up2.fq.gz --failed_out fout.fq.gz

```

Ejemplos 
```python
!fastp -i /content/drive/MyDrive/PR69/HA1AB3SS04_S4_L1_R1_001.fastq.gz -I /content/drive/MyDrive/PR69/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o /content/drive/MyDrive/PR69/Salidas/FastPR1filtrado.fastq.gz -O /content/drive/MyDrive/PR69/Salidas/FastPR2filtrado.fastq.gz -l 36 --cut_right --cut_front -c -m --merged_out /content/drive/MyDrive/PR69/Salidas/mergedPR69.fastq.gz --unpaired1 /content/drive/MyDrive/PR69/Salidas/Up1PR69.fastq.gz --unpaired2 /content/drive/MyDrive/PR69/Salidas/Up2PR69.fastq.gz --failed_out /content/drive/MyDrive/PR69/Salidas/failedPR69

##Filtrado desduplicado
!fastp -i /content/drive/MyDrive/PR69/Salidas/FastPR1filtrado.fastq.gz -I /content/drive/MyDrive/PR69/Salidas/FastPR2filtrado.fastq.gz -o /content/drive/MyDrive/PR69/Salidas/FastPR1filtradodesduplicado.fastq.gz -O /content/drive/MyDrive/PR69/Salidas/FastPR2filtradodesduplicado.fastq.gz -D --dedup


!fastp -i /content/drive/MyDrive/PR69/Salidas/filtered_tile1_1.fastq -I /content/drive/MyDrive/PR69/Salidas/filtered_tile2_2.fastq -o /content/drive/MyDrive/PR69/Salidas/FastPR1filtrado_bbmap.fastq.gz -O /content/drive/MyDrive/PR69/Salidas/FastPR2filtrado_bbmap.fastq.gz -l 40 --length_limit 280 --cut_right --cut_front -c -y -p -m --merged_out /content/drive/MyDrive/PR69/Salidas/mergedPR69_bbmap.fastq.gz --unpaired1 /content/drive/MyDrive/PR69/Salidas/Up1PR69_bbmap.fastq.gz --unpaired2 /content/drive/MyDrive/PR69/Salidas/Up2PR69_bbmap.fastq.gz --failed_out /content/drive/MyDrive/PR69/Salidas/failedPR69_bbmap

```


# 6.2 Trimmomatic

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

Sin embargo en google colab esto no parece ser necesario puesto que simplemente se llama el shell o bash con la correspondiente función (me falta investigar más a fondo este aspecto). 


```shell
%%shell
trimmomatic PE -phred33 -threads 16 \ fwd \ rev  \ /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001.trim.fq.gz \ /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001.trim.fastq.gz /Analisis_Posdoc/PR69/HA1AB3SS04_S4.trim.fq.gz \ ILLUMINACLIP:${adapters}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:90 CROP:150
```
Es muy importante señalar que en este caso hay que tomar en cuenta el orden en el que se colocan los parámetros.

>ILLUMINACLIP:${adapters} eliminamos adaptadores con cierta frecuencia (2:30:10)

>SLIDINGWINDOW: 4:20 cuatro nucleótidos en promedio tienen una calidad menor a 20 se elimina la secuencia (incluyendo el par).

>MINLEN: mínimo de largo 90

>CROP: Cortar la lectura a una longitud de 150


Trimmomatic realiza un *trimming* de calidad adaptativo, cortado de cabeza y cola y remoción de adaptadores. Se puede revisar la documentación y bajar el programa [aquí](http://www.usadellab.org/cms/index.php?page=trimmomatic).

Una de las ventajas del programa es que permite trabajar con secuencias *Paired end*, reteniendo solamente pares coincidentes.
Otra ventaja es que permite coincidencias parciales y *overlapping* para la búsqueda de adaptadores.

Las opciones que podemos utilizar son las siguientes:

## 6.2.1 Eficiencia y formato

>Las siguientes se usan siempre antes de la invocación de los archivos de entrada y salida

- *threads*: este ajuste modifica el número de "hilos" de CPU que Trimmomatic debería usar en las computaciones. Una computadora típicamente tiene cerca de 2 núcleos. los cuales deberían corresponder a una cantidad de 4 hilos disponibles. 
- *phred*:  [-phred33 	-phred64]:  Este ajuste le dice al programa que codifique el archivo

> A partir de aquí son las opciones que van después de la invocación de *Inputs/Outputs* (estas opciones se presentan en mayúsculas) 

Opciones para cambiar la codificación (ver en [Apartado 5](#5-pre-procesamiento-analisis-de-calidad-usando-fastqc-y-fastp)):
Si se requiere leer la codificación de un tipo y sacar la codificación de uno diferente, éstas opciones son las que se necesitan utilizar.

- TOPHRED33: Convierte *scores* de calidad a Phred-33
- TOPHRED64: Convierte *scores* de calidad a Phred-64 

### 6.2.1.1 Cortado (*Cropping*)

Trimmomatic cuenta con varias opciones que pueden ser usadas simultáneamente o no:

-LEADING: Corta bases del inicio de una lectura, si está por debajo del umbral de calidad - adaptativa
-TRAILING: Corta bases del final de una lectura, si está por debajo del umbral de calidad - adaptativa
-CROP: Corta la lectura a una longitud específica
-HEADCROP: Corta el número específicado de bases del inicio de una lectura.

LEADING y TRAILING son cortado adaptativo, lo que significa que cortarán el inicio/fin de las lecturas si fallan la calidad especificada. Lo anterior difiere de CROP y HEADCROP, los cuales podrían cortar a una longitud o número de bases específicas (respectivamente), en este caso, el programa realizara el corte para todas las lecturas.

-MINLEN: se deshará de todas las lecturas que caen bajo una longitud especificada.

### 6.2.1.2 *Trimming* de calidad adaptativo

-SLIDINGWINDOW: realiza un trimming en una ventana de deslizamiento, cortando una vez que la calidad promedio caiga de un umbral especificado.
Toma dos valores como `SLIDINGWINDOW:4:15` lo que significa "Escanear la lectura con una amplitud de ventana de 4 bases, cortando cuando la calidad promedio por base caiga debajo de 5'   

### 6.2.1.3 *Trimming* de adaptadores

Finalmente, trimmomatic tomará un archivo con las secuencias de los adaptadores y las cortará. Siguiendo por ejemplo la llamada: `ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip treshold>:<simple clip treshold>` dónde:

-fastaWithAdaptersEtc: Especifica el *path* a un archivo fasta conteniendo todos los adaptadores, secuencias PCR etc. El nombre de las diferentes secuencias dentro de este archivo determina como serán usadas.
-seedMismatches: Especifica la cuenta máxima de *mismatches*, lo cuál podría seguir permitiendo una coincidencia completa.
-palindromeClipTreshold: especifica que tan precisa es la coincidencia entre las dos lecturas 'ligadas por adaptador' que deben ser palíndromo para las lecturas *PE*
-simpleClipTreshold: especifica que tan precisa debe ser la coincidencia entre cualquier adaptador, etc contra una lectura.

  ![trimmomatic_adapter](https://user-images.githubusercontent.com/13104654/211375856-b34becba-e0e4-450d-8d0b-06552b13b296.png)

Existen 4 posibles escenarios que Trimmomatic puede cubrir:
A. Una secuencia técnica es completamente cubierta por la lectura y así, un alineamiento simple podrá identificarla.
B. Solamente existe una coincidencia parcial entre la secuencia técnica y la lectura y así, es necesariio un alineamiento corto.
C. y D. Ambos pares son probados a la vez, permitiendo que suceda lo siguiente: "es mucho más confiable que un alineamiento corto (B.) y permite que se detecte la lectura del adaptador incluso cuando solo se ha secuenciado una base de este."

El umbral de clip palindrómico escencialmente dice que tan preciso debe ser el alineamiento de adaptadores. Esto es la probabilidad log10 de obtener una coincidencia por una posibilidad aleatoria, y así, los valores alrededor de 30 son recomendados.

Referencias: https://jshleap.github.io/bioinformatics/writting-jNGS_tutorial/#encoding
[Bolger *et al.*, 2014](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)

# 6.3 Otras herramientas
## 6.3.1 Cutadapt

Cutadapt busca el adaptador en todas las lecturas y lo remueve cuando lo encuentra. A menos que se use la opción de filtrado, todas las lecturas que están presentes en el archivo *input* estarán presentes en el *output*, algunas de ellas ya con trimm, otras no. Incluso las lecturas que fueron coradas a una longitud de 0 son un output. Esto puede ser modificado en el comando (opciones).
Puede detectar múltiples tipos de adaptadores. Adaptadores 5' preceden la secuencia de interés mientras que los 3' la siguen. Las distinciones se hacen dependiendo de donde se presenta la secuencia en la lectura. Además también permite el procesamiento de lecturas *Paired End*.


![imagen](https://user-images.githubusercontent.com/13104654/212785938-7cfd92e4-acfd-4e47-86af-cc90123776c3.png)


Con lo anterior en mente, es una herramienta versátil para remover *primers* o en general oligos de las regiones que flanquean el DNA. La guía de uso se puede encontrar en su [página](https://cutadapt.readthedocs.io/en/stable/guide.html).

con el siguiente bloque de código se puede utilizar

```bash

cutadapt -a ADAPTER_FWD -A ADAPTER_REV -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

```
Referencia: https://jshleap.github.io/bioinformatics/writting-jNGS_tutorial/#encoding

## 6.3.2 Seqkit

[Seqkit](https://bioinf.shenwei.me/seqkit/) es un programa o línea de comandos, que permite, no solo eliminar duplicados, si no también manipular secuencias (solamente en formato fastq) eficientemente. El código fuente se puede encontrar [aquí](https://github.com/shenwei356/seqkit).  

> (formatos BAM y SAM tienen que ser convertidos a fastq antes de utilizarse con seqkit).

El siguiente bloque de código funciona localmente para eliminar duplicados

```Bash
seqkit rmdup -s -o clean.fastq input.fastq 
#Remueve secuencias duplicadas del input.fastq y las guarda en clean.fastq 
```

Posiblemente en google colab
```Bash
!seqkit rmdup -s -o ./sec_limpias.fastq ./archivo_salida_dePreproc/PR69.fastq
```
> Se tendría que tomar como input el archivo de salida del último paso del preprocesamiento, que podría ser el de salida de cutadap.

## 6.3.3 Trim-Galore
Para instalar se puede realizar con conda
```bash
conda install -c bioconda trim-galore
```


:alien: 👽 :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien: :alien:

## 6.4 Número y longitud de secuencias después del filtrado de calidad

Una vez realizado el filtrado se pueden correr de nuevo los análisis de calidad (usando fastqc/multiqc o fastp). Además podemos utilizar nuevamente los comandos de bash para analizar longitudes de lecturas, etc.

El bloque de código para estas revisiones (realizando el trimming con trimmomatic) sería el siguiente:

> (También revisamos que se encuentren pareados, esta parte se puede omitir)

El siguiente bloque permite revisar si los archivos están pareados
```bash
##Para explorar que los archivos están pareados
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001_Trim.fastq.gz | wc -l
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001_Trim.fastq.gz | wc -l
```
si queremos saber cual es el número de secuencias usamos el siguiente bloque de código, de nuevo usamos zgrep por que es un archivo comprimido y ahora usaremos los archivos de salida del trimming (o + cutadapt + seqkit).

```bash
##revisar nuevamente la cantidad de secuencias
%%bash
zgrep '^@M02521' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001_Trim.fastq.gz | wc -l
zgrep '^@M02521' /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001_Trim.fastq.gz | wc -l
```

Exploramos la longitud de secuencias. Para ello podemos usar awk de nuevo.

Para cada línea de secuencia podemos contar cada caracter usando el parámetro NR (número de registros) y usando el contador y añadiendo para imprimir en txt.

```bash
%%bash
zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R1_001_Trim.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}}' | sort -n | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/read_length1Trim.txt

zcat /content/drive/MyDrive/Analisis_Posdoc/PR69/HA1AB3SS04_S4_L1_R2_001_Trim.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths){print l, lengths[l]}}' | sort -n | uniq -c > /content/drive/MyDrive/Analisis_Posdoc/read_length2Trim.txt


```

el txt generado lo podemos importar y transformar nuevamente a csv para después graficarlo en matplotlib, se aplican los siguientes bloques:

```python
read_file = pd.read_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_length1Trim.txt',header=None)
read_file.to_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_lengthR1Trim.csv', index=None,header=None)

read_file = pd.read_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_length2Trim.txt',header=None)
read_file.to_csv (r'/content/drive/MyDrive/Analisis_Posdoc/read_lengthR2Trim.csv', index=None,header=None)

```

## 6.5 PhiX 

> Es posible usar esta librería y aplicar este código, sin embargo, en este caso no se hizo una corrida por lo que este código queda en stand by y se podría usar para futuras secuenciaciones (podemos usar el código para evitar remanentes de un control interno de la secuenciación)

Para control de calidad interno 
[Phix](https://support.illumina.com/bulletins/2017/02/what-is-the-phix-control-v3-library-and-what-is-its-function-in-.html) o la librería de control v3 PhiX (FC-110-3001) es derivada del genoma del bacteriófago, bien caracterizado, PhiX. Es una librería concentrada de Illumina (10 nM en 10 µl) que tiene un tamaño promedio de 500 pb y cuenta con una composición balanceada de bases a ~45% GC y ~55% AT.

Puede servir como un control de calibración y puede ser secuenciado solo,  para usarlo en: Cálculos de phasing y prephasing, generación de matrices Cross talk (ref cruzada) y examinar el desempeño promedio de las plataformas de secuenciación Illumina. 

Se puede bajar un programa [BWA](https://github.com/lh3/bwa#type), para realizar el index de la referencia.
BWA: es un paquete de software para el mapeo de secuencias de ADN contra un gran genoma de referencia, como el genoma humano. A mayor número de lecturas de Phix, mayor probabilidad de tener un Profago en el genoma problema. Evidentemente me quedaría con las secuencias que no mapearon vs PhiX.

El siguiente bloque de código es SOLAMENTE para aplicarse en el entorno local de una máquina (no en colab)

```bash
#qc corresponde en este caso a la carpeta de trabajo donde se encuentran todas las secuencias en el presente bloque de código, dependerá de cada  persona y máquina donde se trabaje (es el nombre que cada quien decide)

cd $HOME/Ensamble #Este comando va a corresponder con la carpeta en la que nos encontremos trabajando localmente, en el caso de colab esto no se realiza 
source activate qc #En este comando estamos activando la carpeta qc como la fuente para trabajar
tar -zxvf PhiX_Illumina_RTA.tar.gz && rm PhiX_Illumina_RTA.tar.gz # hay que definir el directorio
bwa index PhiX/Illumina/RTA/Sequence/Chromosomes/phix.fa -p ´01_qc/phix´ # Este es el comando que corresponde al programa BWA para hacer el index
conda deactivate
```
Para colab, primeramente hay que tener considerada la [instalación](#instalacion-de-bwa-y-samtools) y definir el [directorio de las secuencias de la librería Phix](https://drive.google.com/drive/folders/1FqceccL4vIdQRCtta0eGVS_NZiBX55Kz). Podríamos usar el siguiente bloque:

```python
!bwa index /content/drive/MyDrive/Analisis_Posdoc/PhiX/Illumina/RTA/Sequence/Chromosomes/phix.fa -p *./calidad1/phix

```


El siguiente bloque (local) se utiliza para Mapear (por Alineamiento) las *lecturas* versus la referencia, en este caso PhiX:
Consta de tres pasos

1. script usando | Indexar y alinear vs esa referencia
2. Selecciona las secuencias que no mapean vs PhiX
3. Sort=ordenar esas secuencias no mapeadas


```bash 
#De igual forma primero nos movemos a la carpeta ensamble y activamos el entorno del archivo qc
cd $HOME/Ensamble 
source activate qc

#busca el alineamiento, selecciona las secuencias que no mapean y se ordenan dichas secuencias 

bwa mem -t 8 01_qc/phix \
    01_qc/R1.trim.fastq.gz \
    01_qc/R2.trim.fastq.gz | \
    samtools view -bS -f4 - | \
    samtools sort -@ 4 - -o 01_qc/PB69.phix.sorted.bam

samtools fastq \
    -1 01_qc/R1.trim.clean.fastq.gz \
    -2 01_qc/R2.trim.clean.fastq.gz \
    -s 01_qc/S.trim.clean.fastq.gz \
    01_qc/PB69.phix.sorted.bam
```

```bash
ls -ltrh 01_qc/
zcat *./R1.trim.clean.fastq.gz | awk 'END{ print NR/4 }'
zcat *./R2.trim.clean.fastq.gz | awk 'END{ print NR/4 }'
```

El bloque en colab (9min) correspondería a:

```python 

#busca el alineamiento, selecciona las secuencias que no mapean y se ordenan dichas secuencias 
#estableceremos threads en 2 (colab utiliza al menos 2 cores)

!bwa mem -t 2 *./calidad1/phix \
    ./R1.filtrado.fastq.gz \
    ./R2.filtrado.fastq.gz | \
    samtools view -bS -f4 - | \
    samtools sort -@ 2 - -o salidas/PR69.phix.sorted.bam

!samtools fastq \
    -1 salidas/R1.filtradas.clean.fastq.gz \
    -2 salidas/R2.filtradas.clean.fastq.gz \
    -s salidas/S.filtradas.clean.fastq.gz \
    salidas/PR69.phix.sorted.bam
```
> index sirve para indexar la base de datos de la secuencia con el flag -p que define la base de datos de salida

> mem es el algoritmo (Maximal Exact Matches) para el alineamiento contra la referencia, las flags: -t define el número de hilos a usar, luego hay que incluir los archivos del mapeo (fwd y reverse, es importante tomar las secuencias ya filtradas), samtools nos sirve para interactuar con las secuencias, seleccionando y ordenando aquellas que no mapean contra phix y las escribe en el archivo (-o) sorted.bam (cambia el formato), comparando  bam (-b) contra fa (-f), y ordenando (sort) usando 4 hilos (-@). 

> fastq transforma las secuencias bam no mapeadas a fastq en archivos específicos (-1 y -2), así como definir las salidas de lecturas de singletones (-s)

Se puede realizar el Conteo de las secuencias, NUEVAMENTE.
.clean son archivos de salida después de mapear vs PHIX (nosotros establecemos el nombre)

![BAM formato](https://user-images.githubusercontent.com/13104654/216114543-cd5f515e-eb53-4024-9d96-9e2421b34902.png)

```bash
%%bash
zcat R1.trim.clean.fastq.gz | awk 'END{ print NR/4 }'
zcat R2.trim.clean.fastq.gz | awk 'END{ print NR/4 }'
```
Si no se encuentran secuencias de Phix entonces se puede realizar el ensamble con las secuencias filtradas (trimmeadas), de lo contrario se usarán las establecidas en -1 y -2.

https://github.com/Microfred/IntroBioinfo/blob/main/Unidad_3/Readme.md


## 6.6 Determinacion de distribucion de *K-meros*

Aunque hoy en día los ensambladores pueden trabajar distintos tamaños de k-meros, hay que tener en cuenta que elegir un tamaño inapropiado pordría afectar enormemente la calidad de un ensamble. Para ello podemos realizar un escaneo antes, tal como lo vimos en el [apartado 5.2.4 Antes del Filtrado](#524-antes-del-filtrado) de Fastp, este realiza y proporciona un heatmap con un conteo de distribución de 5-meros, lo que en primera instancia nos podría servir para tener una idea de la distribución de *k-meros*. 

(contar k meros, para modificar ciertos parámetros para el ensamble)

Existen otras herramientas que podemos usar  para esta tarea, entre ellas:

### 6.6.1 Kmergenie
[kmergenie](http://kmergenie.bx.psu.edu/) estima la mejor longitud de *K-mero* para en ensamble *De Novo*. Dado un conjunto de lecturas, KmerGenie primero cálcula un histograma de abundancias para muchos valores de *K*. Entonces, para cada valor de *k*, predice el número de distintos *k-meros* genómicos en el dataset y entrega la longitud que maximiza este número. Las precicciones pueden ser aplicadas a ensambles single-k (Velvet, SOAPdenovo 2, ABySS, Minia). Sin embargo ensambladores multi-k se desempeñan mejor, generalmente, con los parámetros por defecto, usando múltiples, más que solo el mejor k predicho. 

El bloque de código de kmergenie en colab, será difícil de aplicar puesto que utiliza instancias de R, lo cual no sé si dificulte la computación. Es posible llamar directamente el repositorio de github para que use esas dependencias. 

> se puede aplicar el ensamble directamente con SPAdes o megahit con los valores por defecto.

### 6.6.2 Velvetadvisor

[Velvet advisor](https://dna.med.monash.edu/~torsten/velvet_advisor/) es una herramienta de cálculo de K-meros, en la cuál, basado en el total de lecturas arroja el valor de K que podría funcionar para el ensamble.
En el caso de PR69 este valor podría ser 287 con 20 veces el *k-mer* coverage para mi ensamble (sugiere entre 10 y 30) 
En el caso de CH606 podría ser de 147 con 20 veces el *k-mer* coverage para mi ensamble (sugiere entre 10 y 30) 

Todos los valores de cobertura en Velvet son proporcionados en cobertura de *k-mer*, por ejemplo que tantas veces tiene que ser visto un *k-mero* entre las lecturas. El índice entre la cobertura de *k-mer* (Ck) y la cobertura estándar (amplitud de nucleótido) (C) es: 

```math
Ck = C * (L - k + 1) / L
```
Dónde `k` es la longitud de hash (pica, trocea, mezcla) y `L` la longitud de la lectura.
La elección de k es un intercambio entre especificidad y sensibilidad. La experiencia muestra que Ck debería encontrarse por encima de 10 para empezar a tener buenos resultados, arriba de 20 se podría estar desperdiciando cobertura. Además, también por experiencia, las pruebas empíricas con diferentes valores de k no son costosas de correr.

![velvet advisor ej](https://user-images.githubusercontent.com/13104654/217328508-8cc231bd-7e81-4355-9312-e773e6278f70.png)

o bien correr las opciones velvetk y [velvet optimizer](https://github.com/tseemann/VelvetOptimiser) para elegir el mejor set.

Además, en caso de utilizar velvet, advisor recomienda usar `-exp_cov auto ` y `-cov_cutoff auto ` en velvetg cuando se exploren los datos por primera vez. 

# 7. Ensamble *De Novo*

Una vez se ha evaluado el control de calidad de las secuencias, éstas se encuentran mezcladas y, como si de un rompecabezas se tratara, hay que armarlo en el orden correcto. El método de alineamiento utilizado para realizar esta tarea se denomina ensamble. Una parte esencial del ensamble es el alineamiento, que involucra disponer de una cantidad masiva de lecturas de DNA, buscando regiones que coincidan unas con otras (regiones de alineamiento) y eventualmente unir el rompecabezas.

Para el Ensamble *De novo*, dependiendo de la plataforma de secuenciación es posible aplicar diferentes estrategias de ensamble:

![4 estrategias de ensamble](https://user-images.githubusercontent.com/13104654/212998923-3620318d-7258-49da-838f-3e63274b195f.png)

[Estrategias de ensamble *De Novo* en secuencias de lectura corta](https://academic.oup.com/bib/article/11/5/457/1746253?login=true)

Los algoritmos de ensamble son la colección de procesos para construir, a partir de cantidades de lecturas de secuencias cortas,  secuencias de DNA original. Las secuencias son alineadas unas con otras y las partes que se superponen son combinadas en una secuencia estrecha. Actualmente, existen dos métodos de algoritmos de ensamble, que diferirán de acuerdo a la complejidad de los datos de secuenciación. 

# 7.1 Métodos
## 7.1.1 OLC (Overlap Layout Consensus) - Ensambladores gráficos de caracteres

Este algoritmo reconoce intersectos entre combinaciones de lecturas para construir una gràfica de las conexiones entre las Lecturas de secuenciación. 
Es una aproximación computacionalmente intensiva donde la complejidad de la computación incrementa con el total de los datos de secuenciación usados durante el ensamblaje. Debido a lo cual este algoritmo se vuelve inexcusable con los secuenciadores como Ilumina, donde millones de lecturas de secuencia corta son necesarias para el ensamble. 

Después de generado el gráfico, visita cada nodo usando un método de la teoría de grafos llamado camino Hamiltoniano, para construir el ensamble final.
La etapa de diseño del algoritmo reduce la complejidad del gráfico preliminar, por condensación de las áreas que surgen sin amigüedad del mismo loci genómico al nodo donde la línea diverge con varios caminos potenciales. 
Es bastante difícil la selección de ese camino. Esto arroja subgrafos para hacer *contigs* que se pueden describir como secuencias *unitig* ensambladas inequívocamente que tienen una alta profundidad de secuenciación y están vinculadas a un gran número de otros *contigs*. Ahora, el unitig se empareja con otros unitigs para formar una secuencia de andamios (*Scaffolds*).

El último paso para realizar el proceso de **consenso** incluye la lectura a través de subgrafos contiguos y extraer la secuencia de consenso para lecturas de cada subgrafo. Otro algoritmo de caracteres involucra la misma teoría de sobrelape de grafos, pero difiere ligeramente ya que simplifica el gráfico removiendo bordes transitivos que tienen detalles redundantes.

Referencias: Li et al., 2012; Chang et al. 2012

![Hamiltonian Path y Consenso](https://user-images.githubusercontent.com/13104654/215958857-f462cfc2-9155-494b-ac6a-13cf6e8d525a.png)
Commins *et al.*, 2009  

![Layout1](https://user-images.githubusercontent.com/13104654/216516419-c8c0b76a-8eae-45cc-ae15-78e1038a0e28.png)

![Layout2](https://user-images.githubusercontent.com/13104654/216521652-b3a7185c-baff-45c0-8531-add7d2a25cb6.png)

![Layout3](https://user-images.githubusercontent.com/13104654/216521784-3320bf35-4111-40f9-9b08-e063bb5363c2.png)

![Consenso](https://user-images.githubusercontent.com/13104654/216521977-9021f8ac-a19a-4ec4-9e85-ff382baa3324.png)


## 7.1.2 Gráficos De Brujin 

Esencialmente, los ensambladores con gráficas De Brujin, rompen las lecturas en subsecuencias de longitud-*k* (*k-meros*), usándolos posteriormente para construir una gráfica. Los nodos de la gráfica representan *k-meros* en este caso, los bordes indican los *k-meros* vecinos que sobrelapan con k-1 base. Las lecturas no se representan como tal si no que se representan por paths. Esta estructura esta basada en la identidad precisa entre los *k-meros*. En la gráfica, los caminos divergentes son generados por errores de secuenciación que reducen la longitud de los *paths*. Los ensambladores usualmente monitorean la cobertura k-mer de cada nodo, lo que permite que el gráfico sea
limpiado eliminando las puntas de baja cobertura. Por lo general, los ensambladores basados en gráficos De Brujin consumen mucha memoria, aunque se utilizan otros métodos para eficientarlo.

![De Brujin](https://user-images.githubusercontent.com/13104654/216419942-5f37e8ee-c7ea-4fa9-a346-a7b9bdaf6811.png)
Ayling *et al.*, 2019

Existen algunas limitaciones inherentes al ensamble con dBg, como la selección inicial del tamaño del *k-mero* con el cuál se construirá la gráfica. Elegir un tamaño inapropiado pordría afectar enormemente la calidad de un ensamble. Pequeños *k-meros* llevan a gráficas más conectadas; los más largos proporcionan mayor especificidad y pocos loops, pero más desconectados como resultado de los gaps o errores dentro de los datos de las lecturas o la falta de cobertura del genoma. Una limitación adicional es que el el tamaño del *k-mero* no puede exceder prácticamente dos menos que la longitud de lectura para generar al menos dos bordes.

![Repeats_ assembly](https://user-images.githubusercontent.com/13104654/216512545-b062af0b-e4e9-4ebd-836f-57ca4ff8974d.png)


![dBGexample](https://user-images.githubusercontent.com/13104654/216523981-95ea24d8-b5a9-4cab-a7f9-1f8256287c0e.png)

![dBGErrorCorrect](https://user-images.githubusercontent.com/13104654/216524134-cacb944c-5fd1-4a42-920c-fea52911c0ce.png)

# 7.2 Ensambladores
> **Resumen de ensambladores en colab**: Los siguientes ensambladores pueden aplicarse en el entorno de colaboratory:  ✅Megahit se aplica sin problemas, ✅SPAdes se aplica sin problemas en colab, ☑️ALGA se realizó solo en WSL, ❌Unicycler no ha funcionado en Colab ni en WSL (probablemente por que usa instancias de SPAdes y no lo puede llamar desde el path de colab, falta revisar como se puede definir el path); ✅Abyss se instaló con conda y se realizó el ensamble sin problema, ✅velvet se pudo instalar con conda pero no se compiló para Kmeros >31 (‼️revisar si es posible compilarlo para mayor número de K-meros), ❌Clover se instaló pero al parecer no funciona. ✅IDBA funcionó sin problema, aunque su desempeño no fue tan bueno. Revisar el ensamble con ‼️Grasshoper, ‼️SOAP de novo y ‼️PlatanusB  en colab. 


## 7.2.1 SPAdes

La información de uso de [SPAdes](https://github.com/ablab/spades) se puede encontrar en el repositorio señalado. Así mismo el artículo que hace referencia a todos los pipelines que involucran el uso de SPAdes es el siguiente: [Prjibelski *et al.*, 2020](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.102)
 
 De manera similar a cualquier otro ensamblador, el objetivo de SPAdes es construir secuencias continuas y precisas (Contigs y Scaffolds) a partir de secuencias cortas. 

SPAdes inicia el pipeline de ensamble construyendo una DBg a partir de las secuencias cortas. Después, la gráfica construida pasa por un procedimiento de simplificación que involucra la eliminación de filos (edges) erróneos. Tales filos son típicamente ocasionados por errores o artefactos de secuenciación. Una vez que la gráfica es simplificada, SPAdes mapea *short paired* y/o *long reads* de nuevo a la gráfica de ensamble usando esta información de alineamiento para realizar resolución repetida o scaffolding usando el módulo exSPAder, el cual construirá paths correctos y continuos para el genoma, siendo ensamblados en la gráfica de ensamble. 
 El protocolo básico 1 es para aislamientos, mientras el protocolo 5 (que puede ser usado más adelante) esta dedicado al descubrimiento de cluster de genes biosintéticos putativos. 
 
 En el caso de colab, primeramente se realiza la instalación pertinente usando conda, tal como en la sección de [instalación de paquetes](#1-instalacion-de-herramientas) en el apartado de [Instalacion de Megahit y SPAdes](#instalacion-de-megahit-y-spades).
 
 La línea de código simplificada para hacer un primer ensamble en colab puede ser:
 
 ```python
 #prueba de instalación
 ! spades.py --test
 
 #abre la ayuda para las flags
 ! spades.py -h
 
 #ensamble de aislamientos
 !spades.py -1 /content/drive/MyDrive/Analisis_Posdoc/PB016/Salidas_Pao/FastPR1_filtrada.fastq.gz -2 /content/drive/MyDrive/Analisis_Posdoc/PB016/Salidas_Pao/FastPR2_filtrada.fastq.gz -t 2 --isolate -o /content/drive/MyDrive/Analisis_Posdoc/PB016/Ensamble_Spades
 ```
 Se usa el flag -isolate para indicar que son aislamientos y que hay bastante profundidad de lo contrario el proceso se detiene (en el caso de la secuenciación de lamgebio si se omite el proceso no se detiene).
 
SPAdes guarda puntos de control para reiniciar desde el último, por lo que para reiniciar se puede utilizar el siguiente bloque, no se debe olvidar colocar la carpeta de salida con la que ya se cuenta para que reanude tomando la información de esta:

```python
!spades.py --restart-from last -t 12 -o /content/drive/MyDrive/Analisis_Posdoc/PB16/Ensamble_SPAdes2 
```
En caso de que se requiera realizar el proceso en dos partes primero se puede usar el flag --only-error-correction, para que haga primero una corrección de errores, luego se puede usar el flag --only-assembler para que realice solo el ensamble y sea menos intensivo el uso de recursos. 

Existen otros flags que se pueden considerar.

## 7.2.2 Megahit
[Megahit](https://github.com/voutcn/megahit#basic-usage)
[Li *et al.*, 2015](https://pubmed.ncbi.nlm.nih.gov/25609793/)

Megahit utiliza gráficos succintos DB (SdBG), los cuales son una representación comprimida de los DBG. Una SdBG codifica una gráfica con m bordes o aristas (edges) un 0(m) bits y soporta 0(1) tiempos transversales de un vórtice a sus vecinos. La implementación de Megahit adiciona un vector-bit de una longitud m para marcar la validez de cada arista (así mismo se respalda la remoción de aristas eficientemente) y un véctor auxiliar de 2kt bits (donde k es el tamaño del k-mero y t es el número de vértices zero-indegree (o aislados)).

 > En teoría de grafos, un vértice o nodo es la unidad fundamental de la que están formados los grafos. Un grafo no dirigido está formado por un conjunto de vértices y un conjunto de aristas (pares no ordenados de vértices), mientras que un grafo dirigido está compuesto por un conjunto de vértices y un conjunto de arcos (pares ordenados de vértices). En este contexto, los vértices son tratados como objetos indivisibles y sin propiedades, aunque puedan tener una estructura adicional dependiendo de la aplicación por la cual se usa el grafo; por ejemplo, una red semántica es un grafo en donde los vértices representan conceptos o clases de objetos. Los dos vértices que conforman una arista se llaman puntos finales ("endpoints", en inglés), y esa arista se dice que es incidente a los vértices. Un vértice w es adyacente a otro vértice v si el grafo contiene una arista (v,w) que los une. La vecindad de un vértice v es un grafo inducido del grafo, formado por todos los vértices adyacentes a v. 

A pesar de las ventajas que estas gráficas representan, no es fácil su construcción por lo que megahit cuenta con un potente algoritmo paralelo para la construcción. Es decir puede explotar el paralelilsmo de las unidades GPU adaptando un algoritmo CX1 BWT

Antes de la construcción de la gráfica, todos los (k+1)-meros de las lecturas input son clasificadas y contadas y solo los (k+1)meros que se presentan al menos d (2 por defecto) veces se mantienen como kmeros sólidos. Este método remueve muchas aristas espurias, pero puede ser riesgoso para el ensamble metagenómico ya que especies de muy baja abundancia pueden haber sido secuenciadas a muy baja profundidad, por lo que también se introduce una estrategia denominada mercy-kmer para recuperar dichos bordes (estos mercy-kmers se agregan a la gráfica para mejorar la contiguidad). 

Se implementa también una estrategia de mútiples tamaños de kmeros, en el cual iterativamente se construyen múltiples SdBGs de un pequeño a un mayor k. Mientras los kmeros pequeños son favorables para filtrar bordes erróneos y rellenar gaps en regiones de baja cobertura, un mayor k es útil para resolver las repeticiones. En cada iteración, se limpian los bordes potencialmente erróneos removiendo puntas, uniendo burbujas y removiendo bordes de baja cobertura local.

 En el caso de colab, primeramente se realiza la instalación pertinente usando conda, tal como en la sección de [instalación de paquetes](#1-instalacion-de-herramientas) en el apartado de [Instalacion de Megahit y SPAdes](#instalacion-de-megahit-y-spades).
 
 
 Básicamente, para el ensamble se usa el siguiente bloque de código para lecturas *Paired end* 
 
 ```Phyton
#Probar lo siguiente para la dirección de las lecturas
#En shell es posible asignar así
#!Lectura1=/content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz
#!Lectura2=/content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz

#!echo $Lectura1
#!echo $Lectura2

#en python podría ser 
lectura1 = path
lectura2 = path

print(lectura1)
print(lectura2)

# !megahit -1 pe_1.fq -2 pe_2.fq -o out #Revisar como se realizará el output en colab
!megahit -1 /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R1_001.fastq.gz -2 /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/HA1AB3SS04_S4_L1_R2_001.fastq.gz -o /content/drive/MyDrive/Analisis_Posdoc/PR69/salidas/ensamble/HA1AB3SS041.megahit_asm

#en caso de que funcione la asignación de variables 
!megahit -1 $Lectura1 -2 $Lectura2 -o ensamble_megahit_res

```
 
 Los contigs pueden encontrarse en el directorio de salida (para colab tenemos que establecerlo) como `final.contigs.fa`

```python
#Crear archivo zip desde carpeta de Colab
!zip -r /content/Ensamble_Megahit.zip /content/Ensamble_Megahit
```

para observar la gráfica de contigs

https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly


```Phyton
!megahit_toolkit contig2fastg 99 k99.contigs.fa > k99.fastg
 ```
 
Algunos tips acerca del ensamble lo podemos encontrar [aquí](https://github.com/voutcn/megahit/wiki/Assembly-Tips)

Para observar el principio del archivo de contigs y contar el número aproximado podemos utilizar el siguiente bloque
```bash

%%bash
head contigs.fa
grep '>' ./directorio_salida/contigs.fa | wc -l 

```
 
## 7.2.3 ALGA 
.
Las herramientas de ensamble basadas en las gráficas de De Brujin son las preferidas para lecturas cortas, debido a que después de la descomposición de las lecturas hay una pérdida de información. Pero también al alto índice de error asociado con nuevas tecnologías que pobremente corresponde con las gráficas de descomposición.
La superioridad de los ensambladores basados en dBG de acuerdo a tiempo y uso de memoria es bien conocido, pero otros ensambladores se desempeñan mejor. Los algoritmos con la estrategia OLC dan contigs más confiables pero con problemas significativos de memoria y tiempo. [ALGA](http://alga.put.poznan.pl/) ha mostrado desempeñarse bien incluso en memoria y tiempo a pesar de ser de tipo OLC ([Swat *et al.*, 2021](https://academic.oup.com/bioinformatics/article/37/12/1644/6104855)).
 
 En el caso de este ensamblador, solo fue posible instalarlo y utilizarlo en WSL:
 
 ```bash 
 #para instalar
wget https://github.com/swacisko/ALGA/archive/refs/tags/1.0.3.tar.gz
tar zvxf 1.0.3.tar.gz

#nos movemos a la carpeta y compilamos
cd 1.0.3
mkdir build
cd build
cmake ..
make

## para ejecutar habría que usar ./ALGA

 ```
### 7.2.3.1 Corrección de errores con Musket
> Esta herramienta recién se ha aplicado en WSL únicamente.

Antes de ensamblar con ALGA se recomienda realizar corrección de errores con musket (se podría usar la corrección de SPAdes?)

Para utilizar [musket](https://musket.sourceforge.net/homepage.html) se recomienda primero su instalación y compilación (ya que esta basado en C) 

```bash
#Instalación
wget https://sourceforge.net/projects/musket/files/musket-1.1.tar.gz/download
tar zvxf musket-1.1.tar.gz

cd musket-1.1

#compilación
mkdir build
cd build
cmake ..
make

#para realizar la corrección de las lecturas
musket -omulti corrected -inorder /mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR1filtradodesduplicado.fastq /mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR2filtradodesduplicado.fastq 

```
### 7.2.3.2 Corrección de errores con SparkEC 
> Aún no se intenta su aplicación

[SparkEC](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05013-1) es una herramienta que trabaja en paralelo, capaz de arreglar aquellos errores generados durante el proceso de secuenciación. Su arquitectura es la misma utilizada por CloudEC (algoritmos MSA), sin embargo se eficientó su desempeño, disminuyendo el tiempo de ejecución y mejorando su utilidad evitando realizar ciertas tareas de forma manual.
Para su instalación y aplicación podemos dirigirnos al [repositorio SparkEC](https://github.com/UDC-GAC/SparkEC)

> Determinar si se puede aplicar el script en colab

Las lecturas ya corregidas se pueden ensamblar con ALGA:

 ```bash
 ./ALGA 
##Primero hay que descomprimir los archivos (ALGA no trabaja con .gz), la descompresión se puede hacer en colab
gzip -d /mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR1filtradodesduplicado.fastq.gz
gzip -d /mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR2filtradodesduplicado.fastq.gz

##luego ensamblamos

./ALGA --file1=/mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR1filtradodesduplicado.fastq --file2=/mnt/c/Users/adria/Downloads/PR69_Ensamble/FastPR2filtradodesduplicado.fastq --threads=14 --output=PR69algacontigs


##Agregar lo siguiente si se piensa que es de baja calidad (probablemente PR69 lo sea) y que 
##aún tiene un alto número de errores

--error-rate=0.02

 ```
 
 
 ## 7.2.4 Unicycler
 
❗Dado que Unicycler utiliza módulos de SPAdes, hasta este momento no se ha logrado (en colab o en WSL) que sea posible ensamblar, quizá se deba al path en el que se encuentra SPAdes 💔 se continuará intentando. 

 ![image](https://user-images.githubusercontent.com/13104654/218617669-7611046b-f3e8-48d2-ad03-46a095c42621.png)
 [Unicycler](https://github.com/rrwick/Unicycler)
 
 
 ## 7.2.5 IDBA
Para instalar en colab:

```python
#Instalar IDBA (ensamblador)
!conda install -c bioconda IDBA -y 
``` 
Antes de comenzar a ensamblar hay que transformar las secuencias de formato fastq a formato fasta 

```python
#Convertir a formato fasta y filtrar N's
!fq2fa --merge --filter /content/drive/MyDrive/PR69/Salidas/MH005_filtered1.fastq /content/drive/MyDrive/PR69/Salidas/MH005_filtered2.fastq /content/drive/MyDrive/PR69/Salidas/MH005_filtered12.fa 

#Ensamble IDBA
!idba_ud -r /content/drive/MyDrive/PR69/Salidas/MH005_filtered12.fa -o /content/drive/MyDrive/PR69/Ensamble_IDBA/_idba_output --mink 51 --maxk 231 --step 10
```


 ## 7.2.6 Velvet
 
 Los primeros bloques permiten la instalación de velvet pero no permiten una modificación para k > 51
 
 ```python
 #Instalar Velvet (ensamblador)
!conda install -c bioconda Velvet -y 
 ```
 
 Para ensamblar en colab 
 ```python
 !velveth /content/drive/MyDrive/PR69/PR69_velvet 51 -fastq -separate -shortPaired /content/drive/MyDrive/PR69/Salidas/FastPR1filtradodesduplicado.fastq.gz /content/drive/MyDrive/PR69/Salidas/FastPR2filtradodesduplicado.fastq.gz
 ```
 
 ❗❗ Buscar la manera de instalar y compilar modificando k 
 
 ## 7.2.7 PlatanusB
 https://github.com/rkajitani/Platanus_B
 
 ## 7.2.8 MeDuSa
 
 
 ## 7.2.9 Abyss
 
 [Repositorio Abyss](https://github.com/bcgsc/abyss#install-abyss-using-conda-recommended)
 
 ```bash
! conda install -c bioconda abyss
 
 # para ambientes dedicados 
! conda create -n abyss-env
! conda activate abyss-env
! conda install -c bioconda abyss
 
 export TMPDIR=/var/tmp
 
 abyss-pe k=25 name=test B=1G \
	in='test-data/reads1.fastq test-data/reads2.fastq'
  
  
  abyss-pe name=ecoli k=96 B=2G in='reads1.fa reads2.fa'


abyss-pe k=96 B=2G name=ecoli lib='pea peb' mp='mpc mpd' \
	pea='pea_1.fa pea_2.fa' peb='peb_1.fa peb_2.fa' \
	mpc='mpc_1.fa mpc_2.fa' mpd='mpd_1.fa mpd_2.fa'

 ```
   Long-distance mate-pair libraries may be used to scaffold an assembly. Specify the names of the mate-pair libraries using the parameter mp. The scaffolds will be stored in the file ${name}-scaffolds.fa. Here's an example of assembling a data set with two paired-end libraries and two mate-pair libraries. Note that the names of the libraries (pea, peb, mpa, mpb) are arbitrary. 
   
 ## 7.2.10 GrassHopper
 
 [Grasshopper](https://sourceforge.net/projects/grasshopper-assembler/)
 [Swiercz *et al.*, 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0202355)
 
 # 7.3 Calidad de Ensamble

**El problema de la validación del genoma**

Aunque se sabe que es muy posible que haya errores en un ensamble de genoma, existe una falta de programas que detecten automáticamente dichos errores o que asignen un *score* de confianza en diferentes regiones del mismo. 

Este problema es particularmente importante ya que la aplicación de las tecnologías de secuencia paralelas, casi siempre generan lecturas cortas y así, incrementan el riesgo de un ensamblado incorrecto (*misassembly*).

Un gran número de estos errores son causados por repeticiones. Los ensambladores pueden ser confundidos por pseudo-sobrelapes entre estas lecturas repetidas (de diferentes copias de repeticiones casi idénticas) y colocarlas juntas. Típicamente pueden inducirse dos tipos de *misassemblies*: el colapso de repeticiones y el rearreglo a gran escala.

La identificación y separación de estas repeticiones colapsadas ha sido estudiada como un problema computacionalmente independiente **el problema de separación de repeticiones** y se han propuesto diversas estrategias combinatoriales y probabilísticas para resolverlo, además el uso de lecturas pareadas puede mejorar esta separación.

El rompecabezas del ensamble casi siempre contiene muchas piezas que son similares en color y forma (repeticiones) y sin tener una idea de como se verá al final. Los errores de ensamble pueden pensarse como piezas que son forzadas a unirse pero que no encajan al final. Se puede definir el termino encajar en el sentido categórico y probabilístico, así, en el caso categórico, los errores de ensamble pueden ser identificados por secuencias que no pueden ser colocadas en el genoma, estas representan secuencias singletones, pares cuya colocación es inconsistente con la librería o sobrelapes cuya composición difiere más de lo que puede ser explicado por errores de secuenciación. En el sentido probabilístico, los errores de ensamble corresponden a regiones del genoma donde el tejido del shotgun es inconsistente con el proceso aleatorio usado para generar  ese secuenciado. Por ejemplo las secciones de un ensamble donde las lecturas se “amontonan” más de lo esperado puede indicar el colapso (coensamblaje) de múltiples copias de una repetición genómica. Este ajuste probabilístico conduce a un elegante formulación del ensamble del genoma como la tarea de identificar un mosaico de lecturas que mejor coincidan con las propiedades del proceso aleatorio utilizado para generar los datos.

El aseguramiento de la calidad aún es más difícil debido al avance continuo de las herramientas de secuenciación, así que cada vez que los datos cambian, representa un nuevo problema para la programación

En ausencia de un genoma de referencia de alta calidad, los nuevos ensambles casi siempre se evaluan dependiendo del número de scaffolds y contigs requeridos para representar el genoma, la proporción de lecturas que puede ser ensamblada, la longitud absoluta de los contigs y scaffolds y la longitud de los contigs y scaffolds relativos al tamaño del genoma.

## Criterios 3C
### Contiguidad
La métrica más comúnmente utilizada es el 💡 **N50**, si todos los contigs en un ensamble se ordenan por longitud, esta métrica es la longitud del contig  más pequeño al 50% de las bases ensambladas.

![Screenshot 2023-08-07 224846](https://github.com/JannaColt/MineriaGen/assets/13104654/045c49a6-cb37-4aef-9f29-ba80db438387)

🔴⚠️ Sólo indica continuidad de bases.

🔴⚠️ Fácil de manipular, no es una medida de precisión del ensamble, debe usarse con precaución.

🔴⚠️ No es significativa para diferentes tamaños de ensamble (no comparable entre especies incluso en el mismo genoma)

🔴⚠️ Sesgado si se excluyen secuencias cortas (lo que casi siempre ocurre)


Herramientas bien establecidas pueden producir ensambles con ⬆️N50, sin embargo, esto puede alcanzarse removiendo k-meros repetidos de bajo coverage (Sacrificando complejidad por contiguidad). Uno puede extender el N50 pero puede carecer de genes conservados.

**NG50** por otro lado, de forma similar a N50 corresponde al contig más pequeño pero al 50% del tamaño del genoma conocido o estimado, por lo que permite comparaciones significativas entre diferentes ensambles. Se pueden filtrar contigs pequeños sin afectar el valor.

Por su parte **L50** corresponde al conteo de contigs en el 50% del ensamble
 Si graficamos la curva Nx, nos daría una mejor visualización ded la continuidad.

### **C**ompletness
- Revisando el tamaño del ensamble
- Nucleótidos conocidos vs desconocidos (esperamos un ensamble sin N's)
- Genes "núcleo" Aseguramiento cuantitativo del genoma ensamblado basado en expectativas evolutivas informadas de  contenidos de genes casi universaldes de ortólogos de copia única. 
- Contenido de k-meros ensamblados
- Mapeo de lecturas y *coverage* de ensamblado

(podemos usar BUSCO y  Merqury: reference-free quality, completeness, and phasing assessment for genome assemblies)

### Correctness (Fidelidad)
Errores que se presentan en el ensamble, proporción del ensamble que está libre de errores como:

✴️ Indels / SNPs

✴️ Mis-joins

✴️ Compresiones repetidas

✴️ Duplicados innecesarios

✴️ Rearreglos.

Pero se realiza contra una referencia, y a veces no contamos con alguna adecuada. Para ello podemos usar dotplots : 
MUMmer dotplot
Chromeister

Algunos pueden utilizar las siguientes estrategias para validar:
▶️ BUSCO/CEGMA para la búsqueda de los genes núcleo

▶️ Mapear lecturas RNASeq y unigenes derivados del ensamble de trasncriptoma

▶️ Mapear proteínas de especies cercanamente relacionadas    

▶️ Mapear lecturas constituyentes que fueron usadas para formar el ensamble y revisar su profundidad y rastreabilidad

▶️ Distribución de NGx (10, 50, 70, 90, etc)

▶️ Distribución de longitud de contigs

▶️ Revisar la presencia de contigs duplicados y otros contaminantes (la forma más fácil es subir el genoma a NCBI)

▶️ Bases constituyentes del ensamble

[Artículo: *De Novo* Genome assembly: what every biologist should know](http://genetica.uab.cat/makingsensegenomicsdata/MakingSenseGenomicData_Reading.pdf)
 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2397507/pdf/gb-2008-9-3-r55.pdf

 ## 7.3.1 QUAST

```python
#Instalación de Quast
!conda install -c bioconda quast -y

#  Correr Quast para ensambles 
!quast.py -o /content/drive/MyDrive/PR69/Estadistica_Ensamble /content/drive/MyDrive/PR69/Ensamble_SPAdes/contigs.fasta /content/drive/MyDrive/PR69/Ensamble_Megahit4/final.contigs.fa ... <dir de cada ensamble>
```

 ## 7.3.2 CheckM
 

## 7.3.3 BUSCO 

 [BUSCO](https://busco.ezlab.org/) es una herramienta con base en expectativas evolutivamente-informadas del contenido de genes ortólogos single-copy casi-universales, la métrica BUSCO es complementaria a las métricas técnicas como N50.

Primero instalamos con conda tratando de que se instalen todos los programas por defecto que se necesiten (esta instalación no permite que se use busco).
```bash
! conda install -c conda-forge -c bioconda -c defaults busco -y -vv
```

Después hay que hacer lo siguiente.

```bash
%cd /content/drive/MyDrive/Busco 
! ls
```


 ```bash
! git clone https://gitlab.com/ezlab/busco.git
 ```

 ```bash
%cd /content/drive/MyDrive/PR69/Busco/busco
! python3 setup.py install
 ```

 ```bash
pd.__version__ 
 ```

 ```bash
sys.path.append('/content/drive/MyDrive/PR69/Busco')
 ```

 ```bash
%cd /content/drive/MyDrive/PR69/Busco/busco/src
! busco -i /content/drive/MyDrive/PR69/Ensamble_SPAdes_MH005_fd/contigs.fasta -l bacteria_odb10 -o anotbusco_PR69 -m genome
 ```

> [!IMPORTANT]
> Una vez que se ha clonado BUSCO en el drive, y que no se ha borrado, podemos omitir el paso de clonado de repositorio y únicamente hacer los pasos de antes y después.

> [!TIP]
> :bulb: 🔆 Para los que están usando la cuenta del laboratorio, el repositorio clonado de BUSCO se encuentra en el path: `/content/drive/MyDrive/Cynthia/CH619_CC/Anotacion/BUSCO`

Si resulta complicado, se puede realizar la anotación usando el web service de [Galaxy](https://usegalaxy.org/)

### 7.3.3.1 Configuracion en Galaxy

1. En el panel de herramientas (primer panel a la izquierda), buscar **Genomic Analysis** -> **Annotation** -> **BUSCO**
2. En el panel central de parámetros de herramientas subir el ensamble del genoma (este son los contig en formato FASTA) 
   
![image](https://github.com/JannaColt/MineriaGen/assets/13104654/93d5fe5e-96d9-4a49-9000-297c3a0aa1af)

luego en los parámetros seleccionar 

:red_circle: en **Lineage data source** -> **Download lineage data** 

:red_circle: en **Mode** -> **Genome assembly (DNA)**

:red_circle: en **Generate miniprot output** se puede dejar a su elección, si se activa tendremos un output tabular con los genes anotados y su posición (gff)

:red_circle: en  **Use Augustus instead of Metaeuk** seleccionamos Augustus, ya que metaeuk es exclusivo de eucariotas.

:red_circle: en **Auto-detect or select lineage?** seleccionamos auto-detect, si sabemos el linaje podemos colocarlo aquí

:red_circle: en **auto-lineage group'*'** podemos dejar el auto lineage de todos los grupos taxonómicos o el prokaryotes 

:red_circle: en **Which outputs should be generated** hay que señalar todos los output que queremos, los cuales pueden ser: un output con el resumen, una listado con los genes perdidos, un gff con datos de la anotación y una imagen con el resumen de genes núcleo.

3. Una vez establecidos todos los parámetros, podemos correr la herramienta. Si existe algún fallo en algunos de los archivos de salida podemos editar los atributos y presionar auto-detect para que se haga la corrección automática.   

*¿Qué me dice el análisis de busco?*

El programa nos proporciona un aseguramiento de la completitud en términos de contenido degenes esperados de un ensamble o conjunto de genes anotados. Los resultaados son simplificados en categorías de completos y de copia única, completos y duplicados, fragmentados o Busco's perdidos (genes marcadores). 

Los resultados de Busco hacen sentido en el contexto de la biología del organismo. Entendiento que los genes duplicados o perdidos pueden ser de origen técnico o biológico. Por ello un alto nivel de duplicación puede ser explicado por un evento de duplicación reciente (biológicamente hablando) o un ensamble quimérico de haplotipos (técnico).

:high_brightness: **Completos** 
Si encontramos genes completos, ya sea de copia única o duplicados, los BUSCO han coincidido con Score suficiente, dentro del rango de Scores esperados y en longitud en cuanto a los alineamientos del perfil BUSCO. Si un ortólogo no está presente en el input, o está parcialmente presente (altamente fragmentado), y un homólogo de alta identidad está presente en longitud completa, es posible que este homólogo haya sido confundido y erróneamente identificado como el BUSCO completo. Los límites del score están optimizados para minimizar esta posibilidad pero aún puede ocurrir.

:low_brightness: **Fragmentados**

Si encontramos genes fragmentados, las coincidencias BUSCO han *Scoreado* dentro del rango de scores pero no dentro del rango de longitud de alineamientos. En ensambles de genomas esto podría indicar que o el hen está parcialmente presente o que el paso de la búsqueda de la secuencia y predicción de genes ha fallado en producir un modelo de longitud completa del gen incluso pensando que el gen está completamente presente en el ensamble. Algunos otros pueden aún estar completos pero pueden ser muy divergentes o tener estructuras genéticas complejas, haciéndolos muy difícil de localizar y predecir al 100%. 

⁉️ :full_moon: **Missing *(perdidos)***

Esto significa que, o no hubo coincidencias significativas, o que las coincidencias BUSCO se presentaron con un score por debajo del rango de scores para el perfil BUSCO. Esto puede indicar que los ortólogos están perdidos, o que el paso de búsqueda de secuencia falló al identificar cualquier coincidencia significativa, o que el paso de predicción de genes falló al producir incluso un modelo genético parcial que podría haber sido reconocido como una coincidencia BUSCO fragmentada. Algunos missing BUSCOs de los aseguramientos de ensambles de genoma podrían así estar parcialmente presentes, e incluso posible pero difícilmente completos, solo que son demasiado divergentes o tienen muy complejas estructuras genéticas haciendo difícil su localización o predicción correcta o incluso parcial.


### 7.3.3.2 Archivos de salida en Colaboratory

Actualmente el perfil de BUSCOs de colaboratory está limitado, probablemente debido a que se tiene que definir el path hacia Augustus. Ver Nota

![image](https://github.com/JannaColt/MineriaGen/assets/13104654/45d1a0b8-85ed-4757-a95b-e79b8e0ca0ac)


![image](https://github.com/JannaColt/MineriaGen/assets/13104654/c4a0cbc3-5bfe-4cde-9566-97cc85a70832)


### 7.3.3.3 Archivos de salida en Galaxy

Los archivos de salida de Galaxy dependerán de lo que se haya marcado en el último punto del paso 2 en [Configuración en Galaxy](#7331-configuracion-en-galaxy)
 
![image](https://github.com/JannaColt/MineriaGen/assets/13104654/c9b4ee98-c07f-4676-85dc-a15307ffe3a7)


> [!NOTE]
> para definir un path ya que colab no acepta export
en entorno normal de linux-unix
```bash
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/bin:$PATH"
export PATH="/path/to/AUGUSTUS/augustus-3.2.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/path/to/AUGUSTUS/augustus-3.2.3/config/"
```
Probablemente esto sirva:
```bash
! echo $PYTHONPATH
%env PYTHONPATH="$/env/python:/content/gdrive/My Drive/Colab Notebooks/path/to/AUGUSTUS/augustus-3.2.3/config/src"
! echo $PYTHONPATH
```

 
 # 7.4 Metaensamblado
 
 
 ## 7.4.1 MAC
 
 ### MAC
Modelo algebraico de adyacencia y clasificación.
Considerar que hay que tener instalado [Mummer](https://colab.research.google.com/drive/1qHrgEQ-rsSG5IxC_BcHzs15cRPEcwwXs#scrollTo=OdaJWjlE0rho&line=2&uniqifier=1) antes de iniciar 

```python
%cd /content/drive/MyDrive/Analisis_Posdoc/PR69/Ensamble/Metaensamble/MACassembler
! ls
## Clonamos el repositorio con código, en la carpeta MACassembler creada

! git clone https://github.com/bioinfomaticsCSU/MAC.git
# Hay que vincular el directorio donde clonamos el repositorio 
sys.path.append('/content/drive/MyDrive/Analisis_Posdoc/PR69/Ensamble/Metaensamble/MACassembler/MAC')

%cd /content/drive/MyDrive/Analisis_Posdoc/PR69/Ensamble/Metaensamble/MACassembler/MAC
!ls

##compilamos usando g++
%%bash
#g++ MAC2.0.cpp -o MAC2

#o si no funciona 
g++ -std=c++11 -std=gnu++11 MAC2.0.cpp -o MAC2.0

#Para correr el programa: 
#Primero hay que posicionarnos en el directorio que contiene las carpetas input, output y temp (ya que no admite path), (tomar en cuenta que la carpeta
#input debe contener los contigs de referencia y query)
#%cd /content/drive/MyDrive/Analisis_Posdoc/PR69/Ensamble/Metaensamble/MACassembler/MAC
%%bash
./MAC2.0 Mh4finalcontigs.fa SPdcontigs.fa
```
 
 ## 7.4.2 Metassembler
 
 #Primero instalamos las herramientas que usará metassembler (omitir si ya están instaladas)
! conda install -c bioconda mummer -y
! conda install -c bioconda bowtie2 -y
! conda install -c bioconda samtools -y

 # 7.5 Mejoramiento del Ensamble
 
 ## 7.5.1 SASpector

[Artículo](https://academic.oup.com/bioinformatics/article/38/10/2920/6564223)

[Repositorio](https://github.com/LoGT-KULeuven/SASpector)
 
 ## 7.5.1 Gap predict

 [Artículo](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8772386/pdf/nihms-1763064.pdf)

 [Repositorio](https://github.com/bcgsc/GapPredict)
 
 ## 7.5.2 Gap closer

 [Artículo](https://academic.oup.com/gigascience/article/9/9/giaa094/5902284?login=false)

 [Repositorio](https://github.com/BGI-Qingdao/TGS-GapCloser)
 
 ## 7.5.3 Improvement assembly
 
[Artículo](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000083)

[Repositorio](https://github.com/sanger-pathogens/assembly_improvement)
 
 ## 7.5.4 Figbird
 
 [Artículo](https://academic.oup.com/bioinformatics/article/38/15/3717/6613135?login=false)

 [Repositorio](https://github.com/SumitTarafder/Figbird)

## 7.5.5 GFinisher

[artículo](https://www.nature.com/articles/srep34963)

[repositorio](https://gfinisher.sourceforge.net/)

## 7.5.6 Filtrado de contigs
### 7.5.6.1 con BBMAP

```bash
#Instalamos BBMAP
! conda install -c bioconda bbmap -y
```

```Bash
# arriba 1000pb
#Para eliminar los contigs que son menor de 1000pb se puede usar bbmap (hay que instalar bbmap primero)
#en in se coloca la dirección del ensamble que vas a filtrar, en out le vas a poner la dirección donde quieres que quede y el nombre de salida
#en min length es la longitud de los que quedan filtrados

! reformat.sh in=/content/drive/MyDrive/Cynthia/CH230_CC/Ensamble_MEGAHIT/CH230_Ensamble_MEGAHIT_final.contigs.fa out=/content/drive/MyDrive/Cynthia/CH230_CC/Limpieza_Ensambles/Limpieza_Ensamble_MEGAHIT/CH230_filteredMAC1000.fasta minlength=1000

```

### 7.5.6.2 con KBase

[Artículo](https://www.nature.com/articles/nbt.4163)

[repositorio app](https://github.com/kbaseapps/kb_assembly_compare/tree/87b74d0dba8cb19dcbd3e7782d4f73a5eaeff76b/ui/narrative/methods/run_filter_contigs_by_length)

Se necesita registrarse a Kbase
[App](https://kbase.us/applist/apps/kb_assembly_compare/run_filter_contigs_by_length/release)
 
 # 8. ANOTACION
 
 ## 8.1 Prokka
 
 ### 8.1.1 Archivos de salida
 Una vez que Prokka ha terminado de anotar, se examina cada uno de los archivos, lo cual se puede realizar utilizando comandos básicos de python.
 Los archivos de salida son los siguientes:
 
 Los archivos GFF y GBK contienen toda la información acerca de las características anotadas (en diferentes formatos).
 El .txt contiene un resumen del número de características anotadas
 El archivo .faa contiene las secuencias de proteínas de los genes anotados
 El archivo .ffn contiene las secuencias nucleotídicas de los genes anotados

 ```bash
#Instalación de Prokka
# !conda install -c bioconda prokka -y 
!conda install -c conda-forge -c bioconda -c defaults prokka
```
 ```bash
!prokka
```

 ```bash
#Correr Prokka 
!prokka --prefix PR69Prokka --locustag PR69 /content/drive/MyDrive/PR69/Ensamble_SPAdes_MH005_fd/contigs.fasta
``` 
 ```bash
#Crear archivo zip desde carpeta de Colab
!zip -r /content/Prokka_PR69_mejor.zip /content/PR69Prokka
```

### 8.1.2 Visualización de características anotadas usando JBrowse

Una forma de visualizar el draft del genoma es utilizando herramientas como *JBrowse genome viewer*

First, we have to make a JBrowse file. Then, we can view it within Galaxy.

## 8.2 PRODIGAL 

## 8.3 PATRIC

La anotación la podemos realizar en PATRIC o BV-BRC
[Web service](https://www.bv-brc.org/)

## 8.4 RAST

# 9. ANOTACIÓN FUNCIONAL

## 9.1 EggNOG Mapper

[Artículo](https://academic.oup.com/mbe/article/38/12/5825/6379734)

[Servicio online](http://eggnog-mapper.embl.de/)

[Repositorio](https://github.com/eggnogdb/eggnog-mapper)

![Diagrama de *workflow*](https://github.com/JannaColt/MineriaGen/assets/13104654/96932911-9cbe-4a7d-9480-286736727738)

## 9.2 FunMappOne 

[Artículo](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2639-2)

[Repositorio](https://github.com/Greco-Lab/FunMappOne)

## 9.3 MicrobeAnnotator

[Artículo](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03940-5)

[Repositorio](https://github.com/cruizperez/MicrobeAnnotator)

## 9.4 Deep Learning
[Artículo](https://pubmed.ncbi.nlm.nih.gov/37010293/)

![DL-functionalAnnotation](https://github.com/JannaColt/MineriaGen/assets/13104654/6535dca0-8913-44b1-9b54-dbb45b829bc6)

## 9.5 MicroScope

[Artículo](https://academic.oup.com/nar/article/48/D1/D579/5606622)
[Plataforma](https://mage.genoscope.cns.fr/microscope/home/index.php)

# 10. MINERÍA GENÓMICA

## 10.1 AntiSMASH

[Artículo](https://academic.oup.com/nar/article/51/W1/W46/7151336)

[Web Tool](https://antismash.secondarymetabolites.org/#!/start)

## 10.2 DeepBGC

[Artículo](https://academic.oup.com/nar/article/47/18/e110/5545735)

[Repositorio](https://github.com/Merck/deepbgc)

## 10.3 e-DeepBGC

[Artículo](https://www.sciencedirect.com/science/article/abs/pii/S0022283622001772)

 [repositorio]()

## 10.4 DecRippter

[Artículo](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001026)

[Repositorio](https://github.com/Alexamk/decRiPPter)

## 10.5 DeepRipp

[Artículo](https://www.pnas.org/doi/10.1073/pnas.1901493116)

[Repositorio](https://github.com/magarveylab/NLPPrecursor)

## 10.6 GECCO

[Artículo](https://www.biorxiv.org/content/10.1101/2021.05.03.442509v1)

[Repositorio](https://github.com/zellerlab/GECCO)

## 10.7 BigCARP

[Artículo](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011162)

[Repositorio](https://github.com/microsoft/bigcarp)

## 10.8 Otros

[Deep Self supervised Learning](https://www.biorxiv.org/content/10.1101/2022.07.22.500861v1.full.pdf)


# 11. ANÁLISIS FILOGENÉTICOS

## 11.1 PhaME 

[Artículo](https://www.nature.com/articles/s41598-020-58356-1)

[Repositorio](https://github.com/LANL-Bioinformatics/PhaME)

[más...](https://phame.readthedocs.io/en/latest/)

 pipeline: TORMES https://github.com/nmquijada/tormes








Ejemplo Metodología artículo
Experimental Design, Materials and Methods
PR69 was previously isolated and identified with 16SrRNA with GenBank accessionnumber #####. Genomic DNA was extracted using the GenomicDNA Purification Kit(NewEnglandBiolabs.). IlluminaHiSeq4000 paired-end (2×151bp) sequencing of PR69 was performed by Langebio (CINVESTAV). The library was processed using the XXXX Library Preparation Kit (96samples (Illumina,Inc.,SanDiego,CA,USA). Total sequencing reads 5,130,218 of 4,447,316 were mapped. Aftermapping, Sambamba[10] and SAMTools[11] were respectively used to remove duplicated reads and identify variants. 
Thereadswereassembledinto58contigs,aGCcontentof41.60%usingGCEAssembler(version1.2;https://cge.cbs.dtu.dk/services/Assembler/)[12].TheassembleddatawasannotatedusingRASTrapidannotationusingsubsystemtechnologyversion2.0[6–8].BacterialsecondarymetabolitebiosynthesisgeneclusterswereidentifiedandannotatedbyantiSMASHversion5.0and6.0usingassembledfastafileoutput[







