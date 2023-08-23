## PIPELINES PARA LLAMADO DE VARIANTES EN SECUENCIAS ONT

### Configuración previa del cluster

El cluster debe contar con los siguientes módulos:
- singularity
- nextflow
- guppy
- tensorflow

### Configuración de los ambientes de conda
Se debe contar con los siguientes ambientes de conda:
- `conda env create -f conda-env-clair3.yml`
- `conda env create -f conda-env-happy.yml`
- `conda env create -f conda-env-nanopolish.yml`

### Configuración inicial del repositorio

Crear las siguientes carpetas

```
git clone https://github.com/lhlopezm2/TGS-ORIGEN.git
cd TGS-ORIGEN
mkdir alignment
mkdir base-calling
mkdir methylation-calling
mkdir metrics
mkdir raw_reads
mkdir results
mkdir variant-calling
```

Dependiendo del pipeline que se vaya a ejecutar se deberá crear una carpeta correspondiente dentro de variant-calling y de metrics.

```
cd variant-calling
mkdir vcf-*-*-*
cd ..
cd metrics
mkdir metrics-*-*-*
cd ..
```

Se deben ubicar las lecturas fast5 en la carpeta raw_reads. Para las presentes pruebas preliminar con datos pequeños usaremos el archivo correspondiente al subset 4 de la secuenciación del individuo HG002, Ashkenazim Trio hijo y tiene un peso de aproximadamente 26 GB.

```
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/combined_2018-08-10/combined_2018-08-10.raw_fast5s.4.tar
```

Al descargar este archivo tenemos que descomprimirlo.
```
mkdir raw_reads
tar -xf combined_2018-08-10.raw_fast5s.4.tar -C raw_reads
```

Ahora bien necesitamos el genoma de referencia y el vcf benchmarking para calcular las metricas.

```
wget -O HG002_GRCh38_benchmark.vcf.gz https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz

wget -O HG002_GRCh38_benchmark.vcf.gz.tbi https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

wget -O GRCh38.bed https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

Para utilizar PEPPER-MARGIN-DEEPVARIANT se debe descargar la imagen de singularity.

```
singularity pull docker://kishwars/pepper_deepvariant:r0.8
```

Para utilizar wf-human-variation se debe clonar el correspondiente repositorio
```
git clone https://github.com/epi2me-labs/wf-human-variation.git
```

Ahora bien, si los modulos y ambientes de conda coinciden con los nombrados al interior de los pipelines, se puede proceder a correr cualquiera de ellos usando slurm.

```
sbatch setup-*-*-*.sh
```

### Explicación de los principales parámetros en las diferentas herramientas bioinformáticas


