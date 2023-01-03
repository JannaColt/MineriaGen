# MINERÍA GENÓMICA 
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

### b) Revisión de las secuencias, número, longitud, que ambos estén pareados (Esta parte se puede omitir ya que fastqc te muestra algunos de estos resultados)

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

![P69 calidad de secuencia por flowcell](https://user-images.githubusercontent.com/13104654/210288959-a6367307-2827-4ff1-b6ec-2dfb468564c0.png)

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


En este apartado se puede generar una  &#x1F4D9;**alerta**, si el contenido de bases varía más del 10% en cualquier posición, y generará un &#x1F534; **fail** si este porcentaje de variación es mayor al 20%.

La causa de este sesgo puede ser el paso de *priming* aleatorio en la producción de las librerías. A pesar de que los hexameros en el priming deben presentarse con igual frecuencia en el mix y deberían realizar el *prime* con eficiencia similar, en la realidad no se da el caso y ciertos hexámeros son favorecidos durante este paso.

¿Entonces, el sesgo tendría implicaciones en los análisis downstream?. Hay algunos puntos a tomar en cuenta:

- Es posible que como parte del sesgo haya un incremento en el *mis-priming* - ocasionando un número alto de *mis-called* bases al inicio de la secuencia, y, 
- Es posible que la selección del sesgo introducido tenga un efecto significativo en la capacidad de la librería de medir el contenido original debido a ciertas secuencias favorecidas.
Sin embargo estos puntos pueden no representar un gran problema ya que son fácilmente detectados,  algunos mencionan que pueden mitigarse por un *Trimming 5'*, sin embargo esto no es un arreglo. Ya que la composición sesgada es creada por la selección de fragmentos de secuenciado y no por errores de llamadas de bases, el único efecto del *trimming* es cambiar de tener una libreria que inicia en posiciones sesgadas a una que inicia más allá de dichas posiciones. 

La única forma de resolver este problema sería introducir nuevos kits de preparación de librerías con una menor disposición al sesgo en el paso del *priming*, sin embargo, a pesar de la advertencia,  no parece que haya consecuencias serias para los análisis posteriores, irónicamente en RNA-seq son más sospechosas las librerías que no presentan este artefacto.

En este [artículo](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085583) se documenta el *mis-priming* en RNA-seq

Algunas de las razones más comunes de este error son:
Secuencias sobrerepresentadas, sesgo en la fragmentación y composición sesgada de las librerías (que a veces ocurre naturalmente).
En el presente caso también se observa una desviación al final, si se está analizando una biblioteca que ha sido recortada agresivamente por el adaptador, naturalmente introducirá un sesgo de composición al final de las lecturas, ya que las secuencias que coinciden con tramos cortos del adaptador se eliminan, dejando solo las secuencias que no coinciden. Por lo tanto, es probable que las desviaciones repentinas en la composición al final de las bibliotecas que han sufrido un recorte agresivo sean falsas.


## 5.1.6.  Contenido de GC por secuencia

![Captura de pantalla 2023-01-02 a la(s) 20 09 24](https://user-images.githubusercontent.com/13104654/210295287-8ed4d1b0-051f-498d-b5f9-d2cdbb5d60c3.png)



## 5.1.7.  Contenido de N por base 

![Captura de pantalla 2023-01-02 a la(s) 20 09 41](https://user-images.githubusercontent.com/13104654/210295307-9831c8c6-f696-4640-ad40-0efc91a27528.png)


## 5.1.8.  Distribución de la longitud de secuencia

![Captura de pantalla 2023-01-02 a la(s) 20 09 57](https://user-images.githubusercontent.com/13104654/210295322-a0759d95-e0a6-42cc-b699-343741a5974e.png)


## 5.1.9.  Niveles de duplicación de secuencias

![Captura de pantalla 2023-01-02 a la(s) 20 10 17](https://user-images.githubusercontent.com/13104654/210295333-f2d9ca7d-0193-4483-b601-e077f178b6c3.png)


## 5.1.10. Secuencias sobre representadas




## 5.1.11. Contenido de adaptadores 

![Captura de pantalla 2023-01-02 a la(s) 20 10 33](https://user-images.githubusercontent.com/13104654/210295352-5c134059-dc4a-4e72-bee1-18e3ee3eadbb.png)





