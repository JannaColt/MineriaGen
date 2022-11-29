# MineriaGen
Contiene pipeline para ensamble y minería genómica de shotgun sequence con google colab

Desde el cuaderno establecido primeramente se instalan todos los paquetes que se usarán y al final se monta el drive en el que se estará trabajando. Es preferible que esto se haga desde el inicio ya que cuando se instala un nuevo paquete se reinicia el entorno y lo que anteriormente llamamos ya no estará disponible. 


Instalamos conda y lo llamamos para proceder con la instalación de los demás paquetes usando conda
```
!pip install -q condacolab
import condacolab
condacolab.install()
```
Luego instalamos miniconda, para lo cual utilizamos un bloque que llame al shell: %%Shell, con lo cual estaremos utilizando los comandos que se usan cuando se utiliza el shell

```
%%shell
#siguiendo el pipeline de FranciscoZorrilla Metagem, el cual es para establecer relaciones metabolicas en metagenomas,
#sin embargo, los primeros pasos para el filtrado de calidad se pueden seguir tal cual, identificar que es posible 
# en el caso de WGS. Por lo pronto habria que instalar paquetes como fastqc y trimmomatic


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

Para montar el drive utilizamos el siguiente bloque:
```
#Montar el drive para trabajar con los archivos dentro del drive
from google.colab import drive
drive.mount('/content/drive')

```
