#!/bin/bash
# Название расчитываемой задачи. Может быть любым.
#SBATCH --job-name="sic"
#
# Множество вычислительных узлов для расчета задачи. Определяет характеристику
# вычислительных узлов.
#SBATCH --partition=intelv4-batch
#
# Запускать каждый расчет на одном узле.
#SBATCH --nodes=1
#
# Расчетное время, после истечения которого задача будет принудительно
# остановлена. В данном случае --- 7 дней.
#SBATCH --time=31-00:00:00
#
# Количество потоков одного процессора (20 для intelv3-batch, 24 для
# intelv4-batch, 256 для knl-batch).
#SBATCH --ntasks-per-node=4

# Чтобы srun заработал с impi
export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

material_name=$1
material=$2

/opt/miniconda3/bin/python main.py $material_name $material
