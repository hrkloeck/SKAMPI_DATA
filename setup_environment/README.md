# SKAMPI_DATA setting up the environment

The example uses [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#) as a work environment and we assue that's already installed on your machine.

Please follow the next steps to setup the environment  

```
git clone https://github.com/hrkloeck/SKAMPI_DATA.git

cd SKAMPI_DATA/setup_environment

conda create --name SKAMPIDATA python=3.9 pip

conda activate SKAMPIDATA

pip install -r requirements.txt_082022
```

Once that is done you may want to run the first example script.
